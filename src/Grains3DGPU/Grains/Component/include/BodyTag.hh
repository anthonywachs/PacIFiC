#ifndef _BODYTAG_HH_
#define _BODYTAG_HH_

#include "Basic.hh"
#include "BitPacker.hh"

// =================================================================================================
/** @brief Body tag bit layout and helper functions for composite-body encoding.

    A body tag is a single uint32 that encodes three fields:
      - Field 0: shapeId (10 bits)          -- index into the RigidBody prototype array (0--1023)
      - Field 1: compositeIdx (14 bits)     -- 0 = standalone; 1--16383 = owning composite index
      - Field 2: subBodyLocalIdx (8 bits)   -- 0 = master / composite CoM; 1--255 = sub-body slot

    A standalone particle has compositeIdx == 0 (and subBodyLocalIdx == 0).
    isSubBody(tag) is equivalent to compositeIdx != 0.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @brief Packer type for the body tag. */
using BodyTagPacker = BitPacker<uint32_t, 10, 14, 8>;

/** @brief Encodes a body tag for a standalone particle.
    @param shapeId Index of the RigidBody prototype this particle uses. */
__HOSTDEVICE__ inline uint makeStandaloneBodyTag(uint shapeId)
{
    BodyTagPacker bp;
    bp.set<0>(shapeId);
    return bp.getValue();
}

/** @brief Encodes a body tag for a composite sub-body.
    @param shapeId       Index of the RigidBody prototype this sub-body uses
    @param compositeIdx  Index of the owning composite (1-based; 0 is reserved for standalone)
    @param localIdx      Local slot within the composite (0 = master) */
__HOSTDEVICE__ inline uint makeSubBodyTag(uint shapeId, uint compositeIdx, uint localIdx)
{
    BodyTagPacker bp;
    bp.set<0>(shapeId);
    bp.set<1>(compositeIdx);
    bp.set<2>(localIdx);
    return bp.getValue();
}

/** @brief Returns the RigidBody prototype index encoded in a body tag. */
__HOSTDEVICE__ inline uint getShapeId(uint bodyTag)
{
    BodyTagPacker bp;
    bp.setValue(bodyTag);
    return bp.get<0>();
}

/** @brief Returns true if the encoded component belongs to a composite. */
__HOSTDEVICE__ inline bool isSubBody(uint bodyTag)
{
    BodyTagPacker bp;
    bp.setValue(bodyTag);
    return bp.get<1>() != 0u;
}

/** @brief Returns the composite index encoded in a sub-body tag (0 = standalone). */
__HOSTDEVICE__ inline uint getCompositeIdx(uint bodyTag)
{
    BodyTagPacker bp;
    bp.setValue(bodyTag);
    return bp.get<1>();
}

/** @brief Returns the sub-body local index encoded in a body tag. */
__HOSTDEVICE__ inline uint getSubBodyLocalIdx(uint bodyTag)
{
    BodyTagPacker bp;
    bp.setValue(bodyTag);
    return bp.get<2>();
}

#endif
