#ifndef _BITPACKER_HH_
#define _BITPACKER_HH_

#include <cstdint>
#include <type_traits>

#include "Basic.hh"

// =================================================================================================
/** @brief The class BitPacker.

    This is a variadic template class that allows packing and unpacking multiple bit fields into a
    single integral type. Each bit field can have a different width specified as template
    parameters.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
template <typename T, size_t... Widths>
class BitPacker
{
    static_assert(std::is_integral_v<T>, "Storage type must be integral");
    static_assert((... + Widths) <= sizeof(T) * 8, "Total widths exceed type size");

private:
    /** @name Parameters */
    //@{
    /** \brief Packed data */
    T data = 0;
    //@}

    /** @name Helpers */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** \brief Helper to get the mask for a specific index */
    template <size_t Index>
    static constexpr T get_mask()
    {
        constexpr size_t width  = get_width<Index>();
        constexpr size_t offset = Offset<Index>::value;
        // Handle the edge case where width is the full size of T to avoid UB
        T mask = (width == sizeof(T) * 8) ? ~T(0) : ((T(1) << width) - 1);
        return mask << offset;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Widths as a constexpr array and compile-time helpers */
    static constexpr size_t widths_arr[] = {Widths...};

    // ---------------------------------------------------------------------------------------------
    /** @brief Helper to get the width for a specific index (constexpr) */
    template <size_t Index>
    static constexpr size_t get_width()
    {
        return widths_arr[Index];
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Nested Offset calculator:
        Offset<0>::value = 0; Offset<N>::value = sum of widths[0..N-1] */
    template <size_t Index, bool IsZero = (Index == 0)>
    struct OffsetImpl
    {
        static constexpr size_t value = OffsetImpl<Index - 1>::value + get_width<Index - 1>();
    };

    template <size_t Index>
    struct OffsetImpl<Index, true>
    {
        static constexpr size_t value = 0;
    };

    template <size_t Index>
    using Offset = OffsetImpl<Index>;
    //@}

public:
    /** @name Constructors and Destructor */
    //@{
    /** @brief Default constructor */
    BitPacker() = default;

    /** @brief Copy constructor */
    BitPacker(const BitPacker&) = default;

    /** @brief Move constructor */
    BitPacker(BitPacker&&) = default;

    /** @brief Copy assignment */
    BitPacker& operator=(const BitPacker&) = default;

    /** @brief Move assignment */
    BitPacker& operator=(BitPacker&&) = default;

    /** @brief Destructor */
    ~BitPacker() = default;
    //@}

    /** @name Get methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Get the entire packed value as raw integer */
    __HOSTDEVICE__ T getValue() const
    {
        return data;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the packed data as raw integer
        @return Quantized integer value */
    template <size_t Index>
    __HOSTDEVICE__ T get() const
    {
        static_assert(Index < sizeof...(Widths), "Index out of bounds for BitPacker::get");
        constexpr T      mask   = get_mask<Index>();
        constexpr size_t offset = Offset<Index>::value;
        return (data & mask) >> offset;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Returns the dequantized floating-point value (fixed-point representation)
        @param min Minimum value in the range
        @param max Maximum value in the range
        @return Dequantized floating-point value */
    template <size_t Index, typename FloatType>
    __HOSTDEVICE__ FloatType getFixed(FloatType min, FloatType max) const
    {
        static_assert(Index < sizeof...(Widths), "Index out of bounds for BitPacker::getFixed");
        constexpr size_t width   = get_width<Index>();
        constexpr T      max_val = (T(1) << width) - 1;

        T quantized = get<Index>();
        if(max <= min)
            max = min + FloatType(1e-10);

        FloatType scale = (max - min) / FloatType(max_val);
        return min + FloatType(quantized) * scale;
    }
    //@}

    /** @name Set methods */
    //@{
    // ---------------------------------------------------------------------------------------------
    /** @brief Set the entire packed value from raw integer
        @param v The complete packed data to set */
    __HOSTDEVICE__
    void setValue(T v)
    {
        data = v;
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Sets a specific part by index with raw integer value
        @param value Raw quantized integer value to set */
    template <size_t Index>
    __HOSTDEVICE__ void set(T value)
    {
        static_assert(Index < sizeof...(Widths), "Index out of bounds for BitPacker::set");
        constexpr T      mask   = get_mask<Index>();
        constexpr size_t offset = Offset<Index>::value;

        // Clear old bits and OR in the new (masked) value
        data = (data & ~mask) | ((value << offset) & mask);
    }

    // ---------------------------------------------------------------------------------------------
    /** @brief Sets a specific part by quantizing a floating-point value (fixed-point)
        @param value Floating-point value to quantize and store
        @param min Minimum value in the range
        @param max Maximum value in the range
        @param saturated Output flag set to true if value was clamped */
    template <size_t Index, typename FloatType>
    __HOSTDEVICE__ void setFixed(FloatType value, FloatType min, FloatType max, bool& saturated)
    {
        static_assert(Index < sizeof...(Widths), "Index out of bounds for BitPacker::setFixed");
        constexpr size_t width   = get_width<Index>();
        constexpr T      max_val = (T(1) << width) - 1;

        saturated = false;
        if(max <= min)
            max = min + FloatType(1e-10);

        FloatType scale = (max - min) / FloatType(max_val);

        if(value <= min)
        {
            set<Index>(0);
            return;
        }

        // Quantize: round to nearest integer
        T quantized = T(llround((value - min) / scale));

        if(quantized > max_val)
        {
            saturated = true;
            quantized = max_val;
        }

        set<Index>(quantized);
    }
    //@}
};

#endif