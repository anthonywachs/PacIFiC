/**
 * @file WriterXML.cpp
 * @brief XML DOM writer utilities built on Xerces-C++ (DOM Level 3 LS).
 *
 * This translation unit provides the `WriterXML` implementation used to create
 * and serialize DOM documents to disk. The writer relies on a shared Xerces
 * runtime lifetime manager (reference-counted) so that multiple XML modules can
 * coexist safely in a single process.
 *
 * The public interface is intentionally minimal: only `WriterXML::m_document`
 * is exposed in the header to preserve legacy expectations.
 */

#include "WriterXML.hh"

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMError.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/dom/DOMException.hpp>
#include <xercesc/dom/DOMLocator.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/util/IOException.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc/util/XMLUniDefs.hpp>

XERCES_CPP_NAMESPACE_USE
using std::string;

/**
 * @brief Current document owned by the writer.
 *
 * The document is created in `initialize()` and released in `terminate()`.
 * The pointer remains `nullptr` when the writer is not initialized.
 */
DOMDocument *WriterXML::m_document = nullptr;

/**
 * @namespace grains_xml_internal
 * @brief Shared Xerces runtime lifetime management.
 *
 * The writer expects these functions to be provided elsewhere (typically in the
 * reader TU) so that ReaderXML and WriterXML can share a single Xerces runtime.
 *
 * @note If the project does not provide shared lifetime management, these can
 *       be replaced by direct `XMLPlatformUtils::Initialize()` and
 *       `XMLPlatformUtils::Terminate()` calls (with appropriate
 * synchronization).
 */
namespace grains_xml_internal {
void xerces_acquire();
void xerces_release();
} // namespace grains_xml_internal

namespace {

/**
 * @brief DOM implementation ID for the DOM Load/Save (“LS”) feature set.
 */
static const XMLCh gLS[] = {chLatin_L, chLatin_S, chNull};

/**
 * @brief Retrieve the Xerces DOM implementation instance for LS.
 *
 * @return DOMImplementation pointer, or null if not available.
 */
DOMImplementation *getImpl() {
  return DOMImplementationRegistry::getDOMImplementation(gLS);
}

/**
 * @brief Downcast a DOM implementation to DOMImplementationLS when supported.
 *
 * @param impl DOM implementation pointer.
 * @return DOMImplementationLS pointer, or null if LS is not supported.
 */
DOMImplementationLS *getImplLS(DOMImplementation *impl) {
  return dynamic_cast<DOMImplementationLS *>(impl);
}

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
   * @brief Transcode a `std::string` to `XMLCh*`.
   * @param s UTF-8 string.
   */
  explicit XStr(const std::string &s) : p(XMLString::transcode(s.c_str())) {}

  /**
   * @brief Transcode a C string to `XMLCh*`.
   * @param s Null-terminated string.
   */
  explicit XStr(const char *s) : p(XMLString::transcode(s)) {}

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
 * @brief Translation-unit local serializer instance.
 *
 * A single serializer is created alongside the document and released at
 * termination. This mirrors legacy behavior where the writer held a global-ish
 * serializer state.
 */
DOMLSSerializer *g_serializer = nullptr;

/**
 * @brief DOM serialization error handler.
 *
 * Xerces may report warnings/errors during serialization. This handler prints
 * diagnostics including URI, line, and column when available.
 *
 * @note Xerces may still abort the write even if `handleError()` returns true.
 */
struct DomErrHandler final : DOMErrorHandler {
  /** @brief Tracks whether any error has been observed. */
  bool had_error = false;

  /**
   * @brief Handle a Xerces DOM serialization error/warning.
   *
   * @param err Error object containing location and message.
   * @return True to request continued processing.
   */
  bool handleError(const DOMError &err) override {
    had_error = true;
    const DOMLocator *loc = err.getLocation();

    std::string uri =
        (loc && loc->getURI()) ? toString(loc->getURI()) : "<unknown>";
    auto line = loc ? loc->getLineNumber() : 0;
    auto col = loc ? loc->getColumnNumber() : 0;

    std::cerr << "Xerces DOM write error"
              << " uri=" << uri << " line=" << line << " col=" << col << " : "
              << toString(err.getMessage()) << "\n";
    return true;
  }

  /**
   * @brief Clear the observed-error state.
   *
   * This function intentionally does not use `override` because
   * `DOMErrorHandler` does not declare `resetErrors()` in all Xerces
   * versions/configurations.
   */
  void resetErrors() { had_error = false; }
};

/** @brief Single TU-local error handler instance attached to the serializer. */
DomErrHandler g_errHandler;

} // namespace

// ----------------------------------------------------------------------------

/**
 * @brief Initialize the XML writer and create a new document.
 *
 * This function:
 * 1. Acquires the shared Xerces runtime,
 * 2. Creates a new DOM document with the provided root element name,
 * 3. Creates a `DOMLSSerializer` and configures pretty-printing and error
 * handling.
 *
 * @param root Name of the document root element.
 * @return Pointer to the document's root element.
 *
 * @pre `root` is not empty (asserted in debug builds).
 * @throws std::runtime_error if the DOM LS implementation is not available.
 *
 * @note The returned root element is owned by the underlying document and
 *       remains valid until `terminate()` releases the document.
 */
DOMElement *WriterXML::initialize(string const &root) {
  assert(root != "");

  grains_xml_internal::xerces_acquire();

  DOMImplementation *impl = getImpl();
  DOMImplementationLS *implLS = getImplLS(impl);
  if (!impl || !implLS) {
    grains_xml_internal::xerces_release();
    throw std::runtime_error("Xerces: DOMImplementation(LS) not available");
  }

  XStr rootName(root);
  m_document = impl->createDocument(nullptr, rootName.get(), nullptr);

  g_serializer = implLS->createLSSerializer();

  // Configure serializer options.
  if (DOMConfiguration *cfg = g_serializer->getDomConfig()) {
    if (cfg->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
      cfg->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);

    // Attach an error handler to surface diagnostics during serialization.
    if (cfg->canSetParameter(XMLUni::fgDOMErrorHandler, &g_errHandler))
      cfg->setParameter(XMLUni::fgDOMErrorHandler, &g_errHandler);
  }

  return m_document->getDocumentElement();
}

/**
 * @brief Serialize the current document to disk and release writer resources.
 *
 * The document and serializer are released regardless of whether serialization
 * succeeds, and the shared Xerces runtime is released after all Xerces-owned
 * objects have been destroyed.
 *
 * @param file Output path for the serialized XML document.
 *
 * @pre `file` is not empty (asserted in debug builds).
 * @pre `initialize()` has been called successfully prior to this function.
 *
 * @note Xerces shutdown must occur after all Xerces objects (serializer,
 * output, document, format targets) have been destroyed. For that reason, the
 * call to `xerces_release()` is performed after the `try` block.
 */
void WriterXML::terminate(std::string const &file) {
  assert(file != "");

  try {
    DOMImplementation *impl = getImpl();
    DOMImplementationLS *implLS = getImplLS(impl);
    if (!impl || !implLS) {
      throw std::runtime_error("Xerces: DOMImplementation(LS) not available");
    }
    if (!m_document || !g_serializer) {
      throw std::runtime_error(
          "WriterXML::terminate called before initialize()");
    }

    // Stream bytes directly to a file target.
    LocalFileFormatTarget target(file.c_str());

    DOMLSOutput *output = implLS->createLSOutput();
    output->setByteStream(&target);

    // Use UTF-8 encoding without allocating transient XMLCh buffers.
    static const XMLCh UTF8[] = {chLatin_U, chLatin_T, chLatin_F, chDigit_8,
                                 chNull};
    output->setEncoding(UTF8);

    // Perform the write.
    g_errHandler.resetErrors();
    g_serializer->write(m_document, output);

    // Release Xerces objects in a safe order.
    output->release();

    g_serializer->release();
    g_serializer = nullptr;

    m_document->release();
    m_document = nullptr;
  } catch (const IOException &e) {
    std::cerr << "Xerces IOException writing " << file << ":\n"
              << "  " << toString(e.getMessage()) << "\n";
  } catch (const XMLException &e) {
    std::cerr << "Xerces XMLException writing " << file << ":\n"
              << "  " << toString(e.getMessage()) << "\n";
  } catch (const DOMException &e) {
    std::cerr << "Xerces DOMException writing " << file << " (code=" << e.code
              << "):\n"
              << "  " << toString(e.msg) << "\n";
  } catch (const std::exception &e) {
    std::cerr << "Exception writing " << file << ":\n"
              << "  " << e.what() << "\n";
  }

  grains_xml_internal::xerces_release();
}

// ----------------------------------------------------------------------------
// Node construction helpers
// ----------------------------------------------------------------------------

/**
 * @brief Create a new element node in the current document.
 *
 * @param name Element tag name.
 * @return Newly created element node.
 *
 * @pre `initialize()` has been called and `m_document` is non-null.
 */
DOMElement *WriterXML::createNode(string const &name) {
  XStr tag(name);
  return m_document->createElement(tag.get());
}

/**
 * @brief Create a new element node and append it to a parent element.
 *
 * @param root Parent element to receive the new child.
 * @param name Child element tag name.
 * @return Newly created child element node.
 *
 * @pre `root` is non-null.
 */
DOMElement *WriterXML::createNode(DOMElement *root, string const &name) {
  DOMElement *node = createNode(name);
  root->appendChild(node);
  return node;
}

/**
 * @brief Set a string attribute on an element.
 *
 * @param root Element whose attribute is set.
 * @param attr Attribute name.
 * @param value Attribute value.
 *
 * @pre `root` is non-null.
 */
void WriterXML::createNodeAttr(DOMElement *root, string const &attr,
                               string const &value) {
  XStr a(attr);
  XStr v(value);
  root->setAttribute(a.get(), v.get());
}

/**
 * @brief Set a numeric attribute on an element.
 *
 * The value is converted using a stream insertion (`operator<<`) into a string
 * and then forwarded to the string-attribute overload.
 *
 * @param root Element whose attribute is set.
 * @param attr Attribute name.
 * @param value Attribute value.
 */
void WriterXML::createNodeAttr(DOMElement *root, string const &attr,
                               double value) {
  std::ostringstream os;
  os << value;
  createNodeAttr(root, attr, os.str());
}

/**
 * @brief Append a text node as the value of an element.
 *
 * @param root Element that receives a text child.
 * @param value Text content.
 *
 * @pre `root` is non-null.
 */
void WriterXML::createNodeValue(DOMElement *root, string const &value) {
  XStr txt(value);
  DOMText *node = m_document->createTextNode(txt.get());
  root->appendChild(node);
}

/**
 * @brief Append a numeric text node as the value of an element.
 *
 * The value is converted using a stream insertion (`operator<<`) into a string
 * and then forwarded to the string-value overload.
 *
 * @param root Element that receives a text child.
 * @param value Numeric content.
 */
void WriterXML::createNodeValue(DOMElement *root, double const &value) {
  std::ostringstream os;
  os << value;
  createNodeValue(root, os.str());
}
