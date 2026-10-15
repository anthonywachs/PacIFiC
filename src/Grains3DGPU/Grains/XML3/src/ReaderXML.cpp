/**
 * @file ReaderXML.cpp
 * @brief XML DOM reader utilities built on Xerces-C++ (DOM Level 3 LS).
 *
 * This translation unit provides the `ReaderXML` implementation and a small
 * amount of process-global Xerces lifetime management.
 *
 * The implementation intentionally preserves legacy behavior where it matters
 * (e.g., node search semantics and certain conversion routines), while using
 * Xerces 3.x APIs (`DOMLSParser`) in place of deprecated types.
 */

#include "ReaderXML.hh"

#include <cassert>
#include <cstdlib> // atoi, atof
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>

#include <xercesc/dom/DOMConfiguration.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/dom/DOMLSParser.hpp>
#include <xercesc/dom/DOMNamedNodeMap.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XMLUniDefs.hpp>

XERCES_CPP_NAMESPACE_USE
using namespace std;

// =====================================================================================
// Shared Xerces init/term (so ReaderXML + WriterXML can both use Xerces safely)
// =====================================================================================

/**
 * @namespace grains_xml_internal
 * @brief Internal Xerces lifetime utilities shared across XML modules.
 *
 * Xerces requires a process-global `Initialize()` / `Terminate()` pairing.
 * This namespace provides a reference-counted acquire/release mechanism to
 * allow multiple translation units (e.g., ReaderXML and WriterXML) to share
 * Xerces safely.
 *
 * @note The symbols in this namespace are intended to remain stable if other
 *       translation units reference them via `extern`.
 */
namespace grains_xml_internal {

/** @brief Guards the refcount and Initialize/Terminate transitions. */
static std::mutex g_mtx;

/** @brief Number of active users of the Xerces runtime. */
static int g_refcount = 0;

/**
 * @brief Acquire the Xerces runtime.
 *
 * Increments a shared reference count and calls
 * `XMLPlatformUtils::Initialize()` on the 0→1 transition.
 *
 * @post Xerces is initialized for the process when this function returns.
 */
void xerces_acquire() {
  std::lock_guard<std::mutex> lock(g_mtx);
  if (g_refcount++ == 0) {
    XMLPlatformUtils::Initialize();
  }
}

/**
 * @brief Release the Xerces runtime.
 *
 * Decrements a shared reference count and calls `XMLPlatformUtils::Terminate()`
 * on the 1→0 transition.
 *
 * @note A non-positive refcount is treated as a no-op to preserve defensive
 *       behavior under mismatched calls.
 */
void xerces_release() {
  std::lock_guard<std::mutex> lock(g_mtx);
  if (g_refcount <= 0)
    return;
  if (--g_refcount == 0) {
    XMLPlatformUtils::Terminate();
  }
}

} // namespace grains_xml_internal

// =====================================================================================
// TU-local state (replaces the old `static DOMBuilder* ReaderXML::m_parser`)
// =====================================================================================
namespace {

/**
 * @brief Translation-unit local DOM parser instance.
 *
 * Xerces exposes `DOMLSParser` as the DOM Level 3 Load/Save parser in
 * Xerces 3.x. A single shared parser is used to match the legacy “global
 * parser” semantics.
 *
 * @warning `DOMLSParser` is not generally documented as thread-safe. This
 * global instance should be treated as single-threaded unless external
 *          synchronization is provided.
 */
DOMLSParser *g_parser = nullptr;

/**
 * @brief DOM implementation ID for the DOM Load/Save (“LS”) feature set.
 *
 * This matches the legacy logic that requested the DOM LS implementation.
 */
static const XMLCh gLS[] = {chLatin_L, chLatin_S, chNull};

/**
 * @brief RAII wrapper for `XMLString::transcode()` allocations.
 *
 * Xerces returns heap-allocated `XMLCh*` buffers from `XMLString::transcode()`.
 * This helper ensures the associated `XMLString::release()` is called.
 */
struct XStr {
  /** @brief Owned Xerces string (released in destructor). */
  XMLCh *p = nullptr;

  /**
   * @brief Transcode an UTF-8 `std::string` to Xerces `XMLCh*`.
   * @param s Input UTF-8 string.
   */
  explicit XStr(const std::string &s) : p(XMLString::transcode(s.c_str())) {}

  /** @brief Release the transcoded buffer. */
  ~XStr() { XMLString::release(&p); }

  /** @brief Access the underlying Xerces string. */
  const XMLCh *get() const { return p; }

  XStr(const XStr &) = delete;
  XStr &operator=(const XStr &) = delete;
};

/**
 * @brief Convert a Xerces `XMLCh*` string to `std::string`.
 *
 * @param x Input Xerces string.
 * @return UTF-8 representation (best-effort) or empty string for null input.
 */
std::string toString(const XMLCh *x) {
  if (!x)
    return {};
  char *c = XMLString::transcode(x);
  std::string s = c ? c : "";
  XMLString::release(&c);
  return s;
}

/**
 * @brief Normalize a DOM node to an element when possible.
 *
 * The legacy code frequently assumes "node-like" pointers are elements.
 * This helper accepts either:
 * - an element node (returned as-is),
 * - a document node (returns its document element),
 * - anything else (returns null).
 *
 * @param n Input DOM node.
 * @return A `DOMElement*` if an element can be produced; otherwise null.
 */
static DOMElement *asElement(DOMNode *n) {
  if (!n)
    return nullptr;

  if (n->getNodeType() == DOMNode::ELEMENT_NODE)
    return static_cast<DOMElement *>(n);

  if (n->getNodeType() == DOMNode::DOCUMENT_NODE) {
    auto *doc = static_cast<DOMDocument *>(n);
    return doc ? doc->getDocumentElement() : nullptr;
  }

  return nullptr;
}

} // namespace

// ----------------------------------------------------------------------------
// Initializes the reader (legacy semantics)
// ----------------------------------------------------------------------------

/**
 * @brief Initialize the XML reader subsystem.
 *
 * This function:
 * 1. Acquires the shared Xerces runtime (reference-counted),
 * 2. Creates the TU-global `DOMLSParser` if it does not already exist,
 * 3. Applies DOM configuration parameters to match legacy behavior.
 *
 * @note Calling `initialize()` multiple times is safe; it is idempotent with
 *       respect to parser creation and increments the Xerces refcount.
 *
 * @throws std::runtime_error if the DOM LS implementation is not available.
 */
void ReaderXML::initialize() {
  grains_xml_internal::xerces_acquire();

  if (!g_parser) {
    DOMImplementation *impl =
        DOMImplementationRegistry::getDOMImplementation(gLS);
    DOMImplementationLS *implLS = dynamic_cast<DOMImplementationLS *>(impl);

    if (!implLS) {
      throw std::runtime_error("Xerces: DOMImplementationLS not available");
    }

    // Xerces 3.x: DOMBuilder is gone; use LSParser instead.
    g_parser =
        implLS->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS, nullptr);

    // Feature configuration mirrors the intent of the legacy DOMBuilder path.
    DOMConfiguration *cfg = g_parser->getDomConfig();

    const bool doNamespaces = false;
    const bool doSchema = false;
    const bool schemaFullChecking = false;

    if (cfg->canSetParameter(XMLUni::fgDOMNamespaces, doNamespaces))
      cfg->setParameter(XMLUni::fgDOMNamespaces, doNamespaces);

    if (cfg->canSetParameter(XMLUni::fgXercesSchema, doSchema))
      cfg->setParameter(XMLUni::fgXercesSchema, doSchema);

    if (cfg->canSetParameter(XMLUni::fgXercesSchemaFullChecking,
                             schemaFullChecking))
      cfg->setParameter(XMLUni::fgXercesSchemaFullChecking, schemaFullChecking);

    // Legacy: fgDOMDatatypeNormalization = true
    if (cfg->canSetParameter(XMLUni::fgDOMDatatypeNormalization, true))
      cfg->setParameter(XMLUni::fgDOMDatatypeNormalization, true);
  }
}

// ----------------------------------------------------------------------------
// Frees the reader (legacy semantics)
// ----------------------------------------------------------------------------

/**
 * @brief Terminate the XML reader subsystem.
 *
 * Releases the TU-global parser (if present) and releases the shared Xerces
 * runtime (reference-counted).
 *
 * @note The legacy code leaked the parser. Releasing it here keeps external
 *       semantics intact while preventing a persistent allocation.
 */
void ReaderXML::terminate() {
  if (g_parser) {
    g_parser->release();
    g_parser = nullptr;
  }

  grains_xml_internal::xerces_release();
}

// ----------------------------------------------------------------------------
// Returns a node from its name and the root node
// ----------------------------------------------------------------------------

/**
 * @brief Return the first descendant element with the given tag name.
 *
 * This matches the legacy “descendant search” behavior by using
 * `getElementsByTagName()` (i.e., not limited to direct children).
 *
 * @param root Root element for the search.
 * @param name Tag name to match.
 * @return First matching node, or null if no match exists.
 */
DOMNode *ReaderXML::getNode(DOMElement *root, string const &name) {
  if (!root)
    return nullptr;
  XStr tag(name);
  DOMNodeList *lst = root->getElementsByTagName(tag.get());
  return lst ? lst->item(0) : nullptr;
}

/**
 * @brief Return the first element child matching a given name.
 *
 * This function iterates over the filtered element-only node list returned by
 * `getNodes(DOMNode*)` and returns the first node whose name matches `name`.
 *
 * @param root Node whose element children are searched.
 * @param name Desired node name.
 * @return Matching node, or null if none exists.
 */
DOMNode *ReaderXML::getNode(DOMNode *root, string const &name) {
  DOMNode *node = nullptr;
  DOMNodeList *nodes = ReaderXML::getNodes(root);

  for (XMLSize_t i = 0; nodes && i < nodes->getLength() && node == nullptr;
       i++) {
    if (name == ReaderXML::getNodeName(nodes->item(i)))
      node = nodes->item(i);
  }
  return node;
}

// ----------------------------------------------------------------------------
// Attribute accessors
// ----------------------------------------------------------------------------

/**
 * @brief Read a scalar attribute and convert it to `double`.
 *
 * Conversion follows legacy semantics via `std::atof()` applied to a transcoded
 * attribute value.
 *
 * @param root Node whose attributes are queried.
 * @param name Attribute name.
 * @return Parsed double value.
 *
 * @pre The attribute exists. Legacy code assumes this and would likely crash.
 *      An `assert()` preserves that expectation in debug builds.
 */
double ReaderXML::getNodeAttr_Double(DOMNode *root, string const &name) {
  DOMNamedNodeMap *nodeValues = root->getAttributes();
  XStr attrName(name);
  DOMNode *value =
      nodeValues ? nodeValues->getNamedItem(attrName.get()) : nullptr;

  assert(value != nullptr);

  const XMLCh *v = value->getNodeValue();
  char *tmp = XMLString::transcode(v);
  double out = std::atof(tmp ? tmp : "");
  XMLString::release(&tmp);
  return out;
}

/**
 * @brief Read an integer attribute and convert it to `int`.
 *
 * Conversion follows legacy semantics via `std::atoi()` applied to a transcoded
 * attribute value.
 *
 * @param root Node whose attributes are queried.
 * @param name Attribute name.
 * @return Parsed integer value.
 *
 * @pre The attribute exists (asserted in debug builds).
 */
int ReaderXML::getNodeAttr_Int(DOMNode *root, const string &name) {
  DOMNamedNodeMap *nodeValues = root->getAttributes();
  XStr attrName(name);
  DOMNode *value =
      nodeValues ? nodeValues->getNamedItem(attrName.get()) : nullptr;

  assert(value != nullptr);

  const XMLCh *v = value->getNodeValue();
  char *tmp = XMLString::transcode(v);
  int out = std::atoi(tmp ? tmp : "");
  XMLString::release(&tmp);
  return out;
}

/**
 * @brief Read an attribute value as a `std::string`.
 *
 * This implementation accepts a `DOMNode*` and requires it to be an element
 * node (or a document node, in which case its document element is used).
 *
 * @param root Node treated as an element.
 * @param name Attribute name.
 * @return Attribute value if present; otherwise an empty string.
 *
 * @note A non-element node is reported to `std::cerr` and yields an empty
 * string, avoiding undefined behavior.
 */
string ReaderXML::getNodeAttr_String(DOMNode *root, string const &name) {
  DOMElement *elem = asElement(root);
  if (!elem) {
    std::cerr << "getNodeAttr_String: node is not an element (type="
              << (root ? root->getNodeType() : -1) << ", attr=" << name
              << ")\n";
    return "";
  }

  DOMNamedNodeMap *attrs = elem->getAttributes();
  if (!attrs)
    return "";

  XStr attrName(name);
  DOMNode *node = attrs->getNamedItem(attrName.get());
  if (!node)
    return "";

  return toString(node->getNodeValue());
}

/**
 * @brief Test whether an attribute exists on a node.
 *
 * @param root Node whose attributes are queried.
 * @param name Attribute name.
 * @return True if the attribute exists; false otherwise.
 */
bool ReaderXML::hasNodeAttr(DOMNode *root, string const &name) {
  DOMNamedNodeMap *nodeValues = root ? root->getAttributes() : nullptr;
  if (!nodeValues)
    return false;

  XStr attrName(name);
  DOMNode *node = nodeValues->getNamedItem(attrName.get());
  return node != nullptr;
}

// ----------------------------------------------------------------------------
// Node name / navigation
// ----------------------------------------------------------------------------

/**
 * @brief Return the node name as a `std::string`.
 *
 * @param root Node pointer.
 * @return Transcoded node name, or empty string for null input.
 */
string ReaderXML::getNodeName(DOMNode const *root) {
  if (!root)
    return "";
  return toString(root->getNodeName());
}

/**
 * @brief Return the “next” node according to legacy semantics.
 *
 * Legacy behavior returns the second entry (`item(1)`) of `getChildNodes()`.
 * This is not a general sibling traversal; it is an index-based selection.
 *
 * @param root Parent node.
 * @return The second child node, or null if missing.
 */
DOMNode *ReaderXML::getNodeNext(DOMNode *root) {
  DOMNodeList *nodes = root ? root->getChildNodes() : nullptr;
  return nodes ? nodes->item(1) : nullptr;
}

// ----------------------------------------------------------------------------
// Node value accessors
// ----------------------------------------------------------------------------

/**
 * @brief Read the text content of a node and convert it to `double`.
 *
 * Legacy behavior reads `root->getFirstChild()->getNodeValue()` and converts it
 * using `std::atof()` after transcoding.
 *
 * @param root Node containing a text child.
 * @return Parsed double value, or 0.0 if the node or its first child is null.
 */
double ReaderXML::getNodeValue_Double(DOMNode *root) {
  if (!root || !root->getFirstChild())
    return 0.0;

  const XMLCh *v = root->getFirstChild()->getNodeValue();
  char *tmp = XMLString::transcode(v);
  double out = std::atof(tmp ? tmp : "");
  XMLString::release(&tmp);
  return out;
}

/**
 * @brief Read the text content of a node as a `std::string`.
 *
 * Legacy behavior reads `root->getFirstChild()->getNodeValue()`.
 *
 * @param root Node containing a text child.
 * @return Text value, or empty string if the node or its first child is null.
 */
string ReaderXML::getNodeValue_String(DOMNode *root) {
  if (!root || !root->getFirstChild())
    return "";
  return toString(root->getFirstChild()->getNodeValue());
}

// ----------------------------------------------------------------------------
// Child/descendant enumeration
// ----------------------------------------------------------------------------

/**
 * @brief Return the list of descendant nodes with the given tag name.
 *
 * This function uses `getElementsByTagName()`, matching legacy semantics of a
 * descendant search rather than direct children only.
 *
 * @param root Root element for the search.
 * @param name Tag name to match.
 * @return Node list of matches, or null if `root` is null.
 */
DOMNodeList *ReaderXML::getNodes(DOMElement *root, string const &name) {
  if (!root)
    return nullptr;
  XStr tag(name);
  return root->getElementsByTagName(tag.get());
}

/**
 * @brief Return a list of element-only children.
 *
 * Legacy behavior clones the input node, removes non-element children from the
 * clone, and returns the clone’s `getChildNodes()`.
 *
 * @param root Node whose children are filtered.
 * @return Node list of element children (backed by a cloned node), or null if
 *         `root` is null.
 *
 * @warning The cloned node is intentionally not released, matching the legacy
 *          behavior that effectively leaked this clone. Callers must assume the
 *          returned `DOMNodeList*` remains valid for the remainder of the
 * process.
 *
 * @note Removing nodes while iterating a live `DOMNodeList` can skip elements.
 *       That quirk is preserved to match historical behavior.
 */
DOMNodeList *ReaderXML::getNodes(DOMNode *root) {
  if (!root)
    return nullptr;

  DOMNode *allNodes = root->cloneNode(true);
  DOMNodeList *nodes = allNodes->getChildNodes();

  for (XMLSize_t i = 0; i < nodes->getLength(); i++) {
    DOMNode *node = nodes->item(i);
    if (node->getNodeType() != DOMNode::ELEMENT_NODE) {
      allNodes->removeChild(node);
    }
  }

  return allNodes->getChildNodes();
}

// ----------------------------------------------------------------------------
// Root node retrieval
// ----------------------------------------------------------------------------

/**
 * @brief Parse an XML document and return its document element.
 *
 * If the shared parser is not initialized, `initialize()` is called implicitly.
 * Parsing uses `DOMLSParser::parseURI()` and returns
 * `doc->getDocumentElement()`.
 *
 * @param xmlFile Path or URI to the XML document.
 * @return Root element of the parsed document, or null on failure.
 *
 * @note Parse errors are reported to `std::cerr` and yield null rather than
 *       throwing to the caller (legacy-friendly behavior).
 */
DOMElement *ReaderXML::getRoot(string const &xmlFile) {
  DOMDocument *doc = nullptr;
  DOMElement *root = nullptr;

  if (!g_parser)
    ReaderXML::initialize();

  try {
    doc = g_parser->parseURI(xmlFile.c_str());
    root = doc ? doc->getDocumentElement() : nullptr;
  } catch (const DOMException &e) {
    XERCES_STD_QUALIFIER cerr << "XML exception " << e.code
                              << XERCES_STD_QUALIFIER endl
                              << toString(e.msg) << XERCES_STD_QUALIFIER endl;
  } catch (const XMLException &e) {
    XERCES_STD_QUALIFIER cerr << "XML exception: " << toString(e.getMessage())
                              << XERCES_STD_QUALIFIER endl;
  }

  return root;
}
