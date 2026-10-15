#ifndef _CONTACTFORCEMODELFACTORY_HH_
#define _CONTACTFORCEMODELFACTORY_HH_

#include "ContactForceModel.hh"
#include "GrainsMemBuffer.hh"

// =================================================================================================
/** @brief The class ContactForceModelFactory.

    Creates the contact force model for each pair of materials.
    This class manages the mapping of material types from strings, as specified in the XML input
    file, to integer IDs, which are then used within the `RigidBody` class.

    1. Material Mapping:
        - Materials are represented as strings in the XML input file, but in the `RigidBody` class,
          materials are stored as integer IDs.
        - To facilitate this, we maintain a map that associates material strings with integer
          material IDs.
        - When constructing a new `RigidBody`, we use this map to ensure that each material string
          is associated with a unique integer ID.
        - If a material string is encountered for the first time, it is added to the map with a new
          ID.
        - The map is stored on the CPU, which works for both CPU and GPU implementations. When
          `RigidBodies` are copied from host to device, the material ID is passed along as a
          parameter.

    2. Handling `contactForceParameters`:
        - The `contactForceParameters` are specified in the XML file using pairs of material
          strings.
        - We store these parameters as a 1D array, where each pair of materials corresponds to an
          index in the array.
        - To access the appropriate `contactForceParameters` for a pair of materials, we use a hash
          function that converts two material IDs into a unique index in the array.
        - For the hash, material IDs for particles and obstacles are considered. We exclude contact
          between two obstacles, but particle-particle and particle-obstacle interactions are
          supported.

    3. Hashing Strategy:
        - We use a triangular indexing approach (via `triangularHash` in GrainsUtils)
          hash = l * ( l + 1 ) / 2 + s, where l is the larger index and s is the smaller one.
          Visualization:
                p0(0)  p1(1)  p2(2)  o0(3)  o1(4)
            p0    0      1      3      6     10
            p1    -      2      4      7     11
            p2    -      -      5      8     12
            o0    -      -      -      9     13
            o1    -      -      -      -     14
          Apparently, indices 9, 13, and 14 will never be accessed, but is fine for now. The better
          way to exclude obst-obst pairs is shown below.
        - Material IDs are assigned such that particle IDs are listed before obstacle IDs. This
          ensures that obstacle-obstacle interactions are not considered in the hash, as they are
          not relevant.
        - The hash is computed by summing the IDs of two materials, with the following table
          serving as an example:
                p0(0)  p1(1)  p2(2)  o0(3)  o1(4)
            p0    0      1      3      6      9
            p1    -      2      4      7     10
            p2    -      -      5      8     11
            o0    -      -      -      -      -
            o1    -      -      -      -      -
          Here, `pX` represents particles, and `oX` represents obstacles. Particle-particle and
          particle-obstacle interactions are supported, while obstacle-obstacle interactions are
          ignored by the hash function. We use a triangular hashing strategy here, and an OFFSET to
          take care of obstacles.

    4. Constraints:
        - **Obstacles Must Be Numbered Last:**
          This is important because we do not consider contact forces between two obstacles. By
          reading particles first and obstacles later, the hash naturally excludes
          obstacle-obstacle interactions.

    5. Error Handling:
        - It is possible for the user to omit certain material pairs from the XML input.
        - We do not validate this during input parsing. If an undefined contact pair is accessed
          during the simulation, it could result in a `NULL` pointer dereference, causing the
          simulation to crash.
        - It is left to the user to ensure that all necessary material pairs are properly defined.

    @author A.YAZDANI - 2024 - Construction */
// =================================================================================================
template <typename T>
class ContactForceModelFactory
{
private:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    ContactForceModelFactory() = default;

    /** @brief Destructor (forbidden) */
    ~ContactForceModelFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    /** @brief Creates and returns the contact force model given an XML node.
        @param root XML node
        @param CF Memory buffer for storing contact force models */
    static void create(DOMElement* root, GrainsMemBuffer<ContactForceModel<T>*, MemType::HOST>& CF);

    /** @brief ContactForceModel objects must be instantiated on device, if we want to use them on
        device. Copying from host is not supported due to runtime polymorphism for this class. This
        function reads a host-side ContactForceModel object, and mimics it in a given device buffer.
        It calls a device kernel that is implemented in the source file.
        @param h_CF Host-side ContactForceModel object
        @param d_CF Device-side ContactForceModel object */
    static void copyHostToDevice(GrainsMemBuffer<ContactForceModel<T>*, MemType::HOST>&   h_CF,
                                 GrainsMemBuffer<ContactForceModel<T>*, MemType::DEVICE>& d_CF);
    //@}
};

#endif
