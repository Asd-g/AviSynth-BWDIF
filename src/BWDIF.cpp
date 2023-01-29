#include <algorithm>

#include "avisynth.h"
#include "BWDIF.h"
#include "VCL2/instrset.h"

template<typename pixel_t, bool spat, int step, int peak, bool debug>
static inline void filterEdge_c(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint_, const int width, const int positiveStride, const int negativeStride, const int stride2,
    float threshold) noexcept
{
    const pixel_t* prev2{ reinterpret_cast<const pixel_t*>(_prev2) };
    const pixel_t* prev{ reinterpret_cast<const pixel_t*>(_prev) };
    const pixel_t* cur{ reinterpret_cast<const pixel_t*>(_cur) };
    const pixel_t* edeint{ reinterpret_cast<const pixel_t*>(edeint_) };
    const pixel_t* next{ reinterpret_cast<const pixel_t*>(_next) };
    const pixel_t* next2{ reinterpret_cast<const pixel_t*>(_next2) };
    pixel_t* __restrict dst{ reinterpret_cast<pixel_t*>(_dst) };

    const pixel_t* prev2Above2{ prev2 - stride2 };
    const pixel_t* prev2Below2{ prev2 + stride2 };
    const pixel_t* prevAbove{ prev + negativeStride };
    const pixel_t* prevBelow{ prev + positiveStride };
    const pixel_t* curAbove{ cur + negativeStride };
    const pixel_t* curBelow{ cur + positiveStride };
    const pixel_t* nextAbove{ next + negativeStride };
    const pixel_t* nextBelow{ next + positiveStride };
    const pixel_t* next2Above2{ next2 - stride2 };
    const pixel_t* next2Below2{ next2 + stride2 };

    typedef typename std::conditional<sizeof(pixel_t) == 4, float, int>::type thresh;
    const thresh thr{ (std::is_integral_v<pixel_t>) ? static_cast<int>(threshold) : threshold };

    for (int x{ 0 }; x < width; ++x)
    {
        if constexpr (std::is_integral_v<pixel_t>)
        {
            const int c{ curAbove[x] };
            const int d{ (prev2[x] + next2[x]) >> 1 };
            const int e{ curBelow[x] };
            const int temporal_diff0{ std::abs(prev2[x] - next2[x]) };
            const int temporal_diff1{ (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) >> 1 };
            const int temporal_diff2{ (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) >> 1 };
            int diff{ std::max({ temporal_diff0 >> 1, temporal_diff1, temporal_diff2 }) };

            if (diff <= thr)
            {
                if constexpr (debug)
                    dst[x] = 0;
                else
                    dst[x] = d;
            }
            else
            {
                if constexpr (debug)
                    dst[x] = peak;
                else
                {
                    if constexpr (spat)
                    {
                        const int b{ ((prev2Above2[x] + next2Above2[x]) >> 1) - c };
                        const int f{ ((prev2Below2[x] + next2Below2[x]) >> 1) - e };
                        const int dc{ d - c };
                        const int de{ d - e };
                        const int maximum{ std::max({ de, dc, std::min(b, f) }) };
                        const int minimum{ std::min({ de, dc, std::max(b, f) }) };
                        diff = std::max({ diff, minimum, -maximum });
                    }

                    dst[x] = std::clamp(std::clamp(!edeint ? (c + e) >> 1 : static_cast<int>(edeint[x]), d - diff, d + diff), 0, peak);
                }
            }
        }
        else
        {
            const float c{ curAbove[x] };
            const float d{ (prev2[x] + next2[x]) * 0.5f };
            const float e{ curBelow[x] };
            const float temporal_diff0{ std::abs(prev2[x] - next2[x]) };
            const float temporal_diff1{ (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) * 0.5f };
            const float temporal_diff2{ (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) * 0.5f };
            float diff{ std::max({ temporal_diff0 * 0.5f, temporal_diff1, temporal_diff2 }) };

            if (diff <= thr)
            {
                if constexpr (debug)
                    dst[x] = 0.0f;
                else
                    dst[x] = d;
            }
            else
            {
                if constexpr (debug)
                    dst[x] = 1.0f;
                else
                {
                    if constexpr (spat)
                    {
                        const float b{ ((prev2Above2[x] + next2Above2[x]) * 0.5f) - c };
                        const float f{ ((prev2Below2[x] + next2Below2[x]) * 0.5f) - e };
                        const float dc{ d - c };
                        const float de{ d - e };
                        const float maximum{ std::max({ de, dc, std::min(b, f) }) };
                        const float minimum{ std::min({ de, dc, std::max(b, f) }) };
                        diff = std::max({ diff, minimum, -maximum });
                    }

                    dst[x] = std::clamp(!edeint ? (c + e) * 0.5f : edeint[x], d - diff, d + diff);
                }
            }
        }
    }
}

template<typename pixel_t, int step, int peak, bool debug>
static inline void filterLine_c(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint_, const int width, const int stride, const int stride2, const int stride3,
    const int stride4, float threshold) noexcept
{
    const pixel_t* prev2{ reinterpret_cast<const pixel_t*>(_prev2) };
    const pixel_t* prev{ reinterpret_cast<const pixel_t*>(_prev) };
    const pixel_t* cur{ reinterpret_cast<const pixel_t*>(_cur) };
    const pixel_t* edeint{ reinterpret_cast<const pixel_t*>(edeint_) };
    const pixel_t* next{ reinterpret_cast<const pixel_t*>(_next) };
    const pixel_t* next2{ reinterpret_cast<const pixel_t*>(_next2) };

    pixel_t* __restrict dst{ reinterpret_cast<pixel_t*>(_dst) };

    const pixel_t* prev2Above4{ prev2 - stride4 };
    const pixel_t* prev2Above2{ prev2 - stride2 };
    const pixel_t* prev2Below2{ prev2 + stride2 };
    const pixel_t* prev2Below4{ prev2 + stride4 };
    const pixel_t* prevAbove{ prev - stride };
    const pixel_t* prevBelow{ prev + stride };
    const pixel_t* curAbove3{ cur - stride3 };
    const pixel_t* curAbove{ cur - stride };
    const pixel_t* curBelow{ cur + stride };
    const pixel_t* curBelow3{ cur + stride3 };
    const pixel_t* nextAbove{ next - stride };
    const pixel_t* nextBelow{ next + stride };
    const pixel_t* next2Above4{ next2 - stride4 };
    const pixel_t* next2Above2{ next2 - stride2 };
    const pixel_t* next2Below2{ next2 + stride2 };
    const pixel_t* next2Below4{ next2 + stride4 };

    typedef typename std::conditional<sizeof(pixel_t) == 4, float, int>::type thresh;
    const thresh thr{ (std::is_integral_v<pixel_t>) ? static_cast<int>(threshold) : threshold };

    for (int x{ 0 }; x < width; ++x)
    {
        if constexpr (std::is_integral_v<pixel_t>)
        {
            const int c{ curAbove[x] };
            const int d{ (prev2[x] + next2[x]) >> 1 };
            const int e{ curBelow[x] };
            const int temporal_diff0{ std::abs(prev2[x] - next2[x]) };
            const int temporal_diff1{ (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) >> 1 };
            const int temporal_diff2{ (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) >> 1 };
            int diff{ std::max({ temporal_diff0 >> 1, temporal_diff1, temporal_diff2 }) };

            if (diff <= thr)
            {
                if constexpr (debug)
                    dst[x] = 0;
                else
                    dst[x] = d;
            }
            else
            {
                if constexpr (debug)
                    dst[x] = peak;
                else
                {
                    const int b{ ((prev2Above2[x] + next2Above2[x]) >> 1) - c };
                    const int f{ ((prev2Below2[x] + next2Below2[x]) >> 1) - e };
                    const int dc{ d - c };
                    const int de{ d - e };
                    const int maximum{ std::max({ de, dc, std::min(b, f) }) };
                    const int minimum{ std::min({ de, dc, std::max(b, f) }) };
                    diff = std::max({ diff, minimum, -maximum });

                    const int interpol{ [&]()
                    {
                        if (!edeint)
                        {
                            if (std::abs(c - e) > temporal_diff0)
                                return std::clamp((((coef_hf[0] * (prev2[x] + next2[x])
                                    - coef_hf[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                                    + coef_hf[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) >> 2)
                                    + coef_lf[0] * (c + e) - coef_lf[1] * (curAbove3[x] + curBelow3[x])) >> 13, d - diff, d + diff);
                            else
                                return std::clamp((coef_sp[0] * (c + e) - coef_sp[1] * (curAbove3[x] + curBelow3[x])) >> 13, d - diff, d + diff);
                        }
                        else
                            return std::clamp(static_cast<int>(edeint[x]), d - diff, d + diff);
                    }() };

                    dst[x] = std::clamp(interpol, 0, peak);
                }
            }
        }
        else
        {
            const float c{ curAbove[x] };
            const float d{ (prev2[x] + next2[x]) * 0.5f };
            const float e{ curBelow[x] };
            const float temporal_diff0{ std::abs(prev2[x] - next2[x]) };
            const float temporal_diff1{ (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) * 0.5f };
            const float temporal_diff2{ (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) * 0.5f };
            float diff{ std::max({ temporal_diff0 * 0.5f, temporal_diff1, temporal_diff2 }) };

            if (diff <= thr)
            {
                if constexpr (debug)
                    dst[x] = 0.0f;
                else
                    dst[x] = d;
            }
            else
            {
                if constexpr (debug)
                    dst[x] = 1.0f;
                else
                {
                    const float b{ ((prev2Above2[x] + next2Above2[x]) * 0.5f) - c };
                    const float f{ ((prev2Below2[x] + next2Below2[x]) * 0.5f) - e };
                    const float dc{ d - c };
                    const float de{ d - e };
                    const float maximum{ std::max({ de, dc, std::min(b, f) }) };
                    const float minimum{ std::min({ de, dc, std::max(b, f) }) };
                    diff = std::max({ diff, minimum, -maximum });

                    const float interpol{ (std::abs(c - e) > temporal_diff0) ? ((coef_hf_f[0] * (prev2[x] + next2[x])
                        - coef_hf_f[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                        + coef_hf_f[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) * 0.25f
                        + coef_lf_f[0] * (c + e) - coef_lf_f[1] * (curAbove3[x] + curBelow3[x])) : (coef_sp_f[0] * (c + e) - coef_sp_f[1] * (curAbove3[x] + curBelow3[x])) };

                    dst[x] = std::clamp(!edeint ? interpol : edeint[x], d - diff, d + diff);
                }
            }
        }
    }
}

/* multiplies and divides a rational number, such as a frame duration, in place and reduces the result */
AVS_FORCEINLINE void muldivRational(int64_t* num, int64_t* den, int64_t mul, int64_t div)
{
    /* do nothing if the rational number is invalid */
    if (!*den)
        return;

    int64_t a;
    int64_t b;
    *num *= mul;
    *den *= div;
    a = *num;
    b = *den;

    while (b != 0)
    {
        int64_t t{ a };
        a = b;
        b = t % b;
    }

    if (a < 0)
        a = -a;

    *num /= a;
    *den /= a;
}

class BWDIF : public GenericVideoFilter
{
    int field_;
    PClip edeint_;
    int opt_;
    float thr_;
    bool pass_;
    bool has_at_least_v8;

    void (*filterEdgeWithSpat)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2,
        float threshold) noexcept;
    void (*filterEdgeWithoutSpat)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2,
        float threshold) noexcept;
    void (*filterLine)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int stride, const int stride2, const int stride3, const int stride4,
        float threshold) noexcept;

    template<typename pixel_t>
    void filter(PVideoFrame& prevFrame, PVideoFrame& curFrame, PVideoFrame& nextFrame, PVideoFrame& dstFrame, PVideoFrame& edeint, const int field, const int tff, IScriptEnvironment* env) noexcept;
    void (BWDIF::* filter_)(PVideoFrame& prevFrame, PVideoFrame& curFrame, PVideoFrame& nextFrame, PVideoFrame& dstFrame, PVideoFrame& edeint, const int field, const int tff, IScriptEnvironment* env) noexcept;

public:
    BWDIF(PClip _child, int field, PClip edeint, int opt, float thr, bool debug, bool pass, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env) override;
    int __stdcall SetCacheHints(int cachehints, int frame_range) override
    {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }
};

template<typename pixel_t>
void BWDIF::filter(PVideoFrame& prevFrame, PVideoFrame& curFrame, PVideoFrame& nextFrame, PVideoFrame& dstFrame, PVideoFrame& edeintFrame, const int field, const int tff, IScriptEnvironment* env) noexcept
{
    int planes_y[3] = { PLANAR_Y, PLANAR_U, PLANAR_V };
    int planes_r[3] = { PLANAR_G, PLANAR_B, PLANAR_R };
    const int* current_planes{ (!vi.IsRGB()) ? planes_y : planes_r };
    const int planecount{ std::min(vi.NumComponents(), 3) };
    for (int i{ 0 }; i < planecount; ++i)
    {
        const int stride{ curFrame->GetPitch(current_planes[i]) / sizeof(pixel_t) };
        const int edeint_stride{ edeintFrame ? edeintFrame->GetPitch(current_planes[i]) / sizeof(pixel_t) : 0 };
        const int dst_stride{ dstFrame->GetPitch(current_planes[i]) / sizeof(pixel_t) };
        const int width{ curFrame->GetRowSize(current_planes[i]) / sizeof(pixel_t) };
        const int height{ curFrame->GetHeight(current_planes[i]) };
        const pixel_t* prev{ reinterpret_cast<const pixel_t*>(prevFrame->GetReadPtr(current_planes[i])) };
        const pixel_t* cur{ reinterpret_cast<const pixel_t*>(curFrame->GetReadPtr(current_planes[i])) };
        const pixel_t* next{ reinterpret_cast<const pixel_t*>(nextFrame->GetReadPtr(current_planes[i])) };
        const pixel_t* edeint{ edeintFrame ? reinterpret_cast<const pixel_t*>(edeintFrame->GetReadPtr(current_planes[i])) : 0 };
        pixel_t* __restrict dst{ reinterpret_cast<pixel_t*>(dstFrame->GetWritePtr(current_planes[i])) };

        env->BitBlt(dstFrame->GetWritePtr(current_planes[i]) + dstFrame->GetPitch(current_planes[i]) * (static_cast<int64_t>(1) - field), dstFrame->GetPitch(current_planes[i]) * 2,
            curFrame->GetReadPtr(current_planes[i]) + curFrame->GetPitch(current_planes[i]) * (static_cast<int64_t>(1) - field), curFrame->GetPitch(current_planes[i]) * 2, curFrame->GetRowSize(current_planes[i]), height / 2);

        prev += static_cast<int64_t>(stride) * field;
        cur += static_cast<int64_t>(stride) * field;
        next += static_cast<int64_t>(stride) * field;
        edeint += static_cast<int64_t>(edeint_stride) * field;
        dst += static_cast<int64_t>(dst_stride) * field;

        const pixel_t* prev2{ (field ^ tff) ? prev : cur };
        const pixel_t* next2{ (field ^ tff) ? cur : next };

        for (int y{ field }; y < height; y += 2)
        {
            if ((y < 4) || (y + 5 > height))
            {
                if ((y < 2) || (y + 3 > height))
                    filterEdgeWithoutSpat(prev2, prev, cur, next, next2, dst, edeint, width,
                        y + 1 < height ? stride : -stride,
                        y > 0 ? -stride : stride,
                        stride * 2, thr_);
                else
                    filterEdgeWithSpat(prev2, prev, cur, next, next2, dst, edeint, width,
                        y + 1 < height ? stride : -stride,
                        y > 0 ? -stride : stride,
                        stride * 2, thr_);
            }
            else
            {
                filterLine(prev2, prev, cur, next, next2, dst, edeint, width,
                    stride, stride * 2, stride * 3, stride * 4, thr_);
            }

            prev2 += stride * static_cast<int64_t>(2);
            prev += stride * static_cast<int64_t>(2);
            cur += stride * static_cast<int64_t>(2);
            next += stride * static_cast<int64_t>(2);
            next2 += stride * static_cast<int64_t>(2);
            edeint += edeint_stride * static_cast<int64_t>(2);
            dst += dst_stride * static_cast<int64_t>(2);
        }
    }
}

BWDIF::BWDIF(PClip _child, int field, PClip edeint, int opt, float thr, bool debug, bool pass, IScriptEnvironment* env)
    : GenericVideoFilter(_child), field_(field), edeint_(edeint), opt_(opt), thr_(thr), pass_(pass)
{
    if (!vi.IsPlanar())
        env->ThrowError("BWDIF: only planar formats are supported.");
    if (vi.height < 4)
        env->ThrowError("BWDIF: height must be greater than or equal to 4.");
    if (field_ < -2 || field_ > 3)
        env->ThrowError("BWDIF: field must be -2, -1, 0, 1, 2, or 3.");
    if (opt_ < -1 || opt_ > 3)
        env->ThrowError("BWDIF: opt must be between -1..3.");

    const int iset{ instrset_detect() };
    if (opt == 1 && iset < 2)
        env->ThrowError("BWDIF: opt=1 requires SSE2.");
    if (opt == 2 && iset < 8)
        env->ThrowError("BWDIF: opt=2 requires AVX2.");
    if (opt == 3 && iset < 10)
        env->ThrowError("BWDIF: opt=3 requires AVX512F.");

    if (thr_ < 0.0f || thr_ > 100.0f)
        env->ThrowError("BWDIF: thr must be between 0.0..100.0.");

    thr_ = (vi.ComponentSize() < 4) ? (thr_ * ((1 << vi.BitsPerComponent()) - 1) / 100.0f + 0.5f) : (thr_ / 100.0f);

    if (edeint_)
    {
        const VideoInfo& vi_ = edeint_->GetVideoInfo();

        if (!vi.IsSameColorspace(vi_))
            env->ThrowError("BWDIF: edeint clip's colorspace doesn't match.");
        if (vi.width != vi_.width || vi.height != vi_.height)
            env->ThrowError("BWDIF: input and edeint must be the same resolution.");
        if (vi.num_frames * ((field_ == -2 || field_ > 1) ? 2 : 1) != vi_.num_frames)
            env->ThrowError("BWDIF: edeint clip's number of frames doesn't match.");
    }

    if ((opt == -1 && iset >= 10) || opt == 3)
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_avx512<uint8_t, true, 32, 255, true> : filterEdge_avx512<uint8_t, true, 32, 255, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_avx512<uint8_t, false, 32, 255, true> : filterEdge_avx512<uint8_t, false, 32, 255, false>;
                filterLine = (debug) ? filterLine_avx512<uint8_t, 16, 255, true> : filterLine_avx512<uint8_t, 16, 255, false>;
                filter_ = &BWDIF::filter<uint8_t>;
            }
            break;
            case 2:
            {
                switch (vi.BitsPerComponent())
                {
                    case 10:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx512<uint16_t, true, 16, 1023, true> : filterEdge_avx512<uint16_t, true, 16, 1023, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx512<uint16_t, false, 16, 1023, true> : filterEdge_avx512<uint16_t, false, 16, 1023, false>;
                        filterLine = (debug) ? filterLine_avx512<uint16_t, 16, 1023, true> : filterLine_avx512<uint16_t, 16, 1023, false>;
                        break;
                    }
                    case 12:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx512<uint16_t, true, 16, 4095, true> : filterEdge_avx512<uint16_t, true, 16, 4095, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx512<uint16_t, false, 16, 4095, true> : filterEdge_avx512<uint16_t, false, 16, 4095, false>;
                        filterLine = (debug) ? filterLine_avx512<uint16_t, 16, 4095, true> : filterLine_avx512<uint16_t, 16, 4095, false>;
                        break;
                    }
                    case 14:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx512<uint16_t, true, 16, 16383, true> : filterEdge_avx512<uint16_t, true, 16, 16383, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx512<uint16_t, false, 16, 16383, true> : filterEdge_avx512<uint16_t, false, 16, 16383, false>;
                        filterLine = (debug) ? filterLine_avx512<uint16_t, 16, 16383, true> : filterLine_avx512<uint16_t, 16, 16383, false>;
                        break;
                    }
                    default:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx512<uint16_t, true, 16, 65535, true> : filterEdge_avx512<uint16_t, true, 16, 65535, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx512<uint16_t, false, 16, 65535, true> : filterEdge_avx512<uint16_t, false, 16, 65535, false>;
                        filterLine = (debug) ? filterLine_avx512<uint16_t, 16, 65535, true> : filterLine_avx512<uint16_t, 16, 65535, false>;
                        break;
                    }
                }

                filter_ = &BWDIF::filter<uint16_t>;
            }
            break;
            default:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_avx512<float, true, 16, 1, true> : filterEdge_avx512<float, true, 16, 1, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_avx512<float, false, 16, 1, true> : filterEdge_avx512<float, false, 16, 1, false>;
                filterLine = (debug) ? filterLine_avx512<float, 16, 1, true> : filterLine_avx512<float, 16, 1, false>;
                filter_ = &BWDIF::filter<float>;
            }
            break;
        }
    }
    else if ((opt == -1 && iset >= 8) || opt == 2)
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_avx2<uint8_t, true, 16, 255, true> : filterEdge_avx2<uint8_t, true, 16, 255, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_avx2<uint8_t, false, 16, 255, true> : filterEdge_avx2<uint8_t, false, 16, 255, false>;
                filterLine = (debug) ? filterLine_avx2<uint8_t, 8, 255, true> : filterLine_avx2<uint8_t, 8, 255, false>;
                filter_ = &BWDIF::filter<uint8_t>;
            }
            break;
            case 2:
            {
                switch (vi.BitsPerComponent())
                {
                    case 10:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx2<uint16_t, true, 8, 1023, true> : filterEdge_avx2<uint16_t, true, 8, 1023, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx2<uint16_t, false, 8, 1023, true> : filterEdge_avx2<uint16_t, false, 8, 1023, false>;
                        filterLine = (debug) ? filterLine_avx2<uint16_t, 8, 1023, true> : filterLine_avx2<uint16_t, 8, 1023, false>;
                        break;
                    }
                    case 12:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx2<uint16_t, true, 8, 4095, true> : filterEdge_avx2<uint16_t, true, 8, 4095, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx2<uint16_t, false, 8, 4095, true> : filterEdge_avx2<uint16_t, false, 8, 4095, false>;
                        filterLine = (debug) ? filterLine_avx2<uint16_t, 8, 4095, true> : filterLine_avx2<uint16_t, 8, 4095, false>;
                        break;
                    }
                    case 14:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx2<uint16_t, true, 8, 16383, true> : filterEdge_avx2<uint16_t, true, 8, 16383, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx2<uint16_t, false, 8, 16383, true> : filterEdge_avx2<uint16_t, false, 8, 16383, false>;
                        filterLine = (debug) ? filterLine_avx2<uint16_t, 8, 16383, true> : filterLine_avx2<uint16_t, 8, 16383, false>;
                        break;
                    }
                    default:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_avx2<uint16_t, true, 8, 65535, true> : filterEdge_avx2<uint16_t, true, 8, 65535, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_avx2<uint16_t, false, 8, 65535, true> : filterEdge_avx2<uint16_t, false, 8, 65535, false>;
                        filterLine = (debug) ? filterLine_avx2<uint16_t, 8, 65535, true> : filterLine_avx2<uint16_t, 8, 65535, false>;
                        break;
                    }
                }

                filter_ = &BWDIF::filter<uint16_t>;
            }
            break;
            default:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_avx2<float, true, 8, 1, true> : filterEdge_avx2<float, true, 8, 1, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_avx2<float, false, 8, 1, true> : filterEdge_avx2<float, false, 8, 1, false>;
                filterLine = (debug) ? filterLine_avx2<float, 8, 1, true> : filterLine_avx2<float, 8, 1, false>;
                filter_ = &BWDIF::filter<float>;
            }
            break;
        }
    }
    else if ((opt == -1 && iset >= 2) || opt == 1)
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_sse2<uint8_t, true, 8, 255, true> : filterEdge_sse2<uint8_t, true, 8, 255, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_sse2<uint8_t, false, 8, 255, true> : filterEdge_sse2<uint8_t, false, 8, 255, false>;
                filterLine = (debug) ? filterLine_sse2<uint8_t, 4, 255, true> : filterLine_sse2<uint8_t, 4, 255, false>;
                filter_ = &BWDIF::filter<uint8_t>;
            }
            break;
            case 2:
            {
                switch (vi.BitsPerComponent())
                {
                    case 10:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_sse2<uint16_t, true, 4, 1023, true> : filterEdge_sse2<uint16_t, true, 4, 1023, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_sse2<uint16_t, false, 4, 1023, true> : filterEdge_sse2<uint16_t, false, 4, 1023, false>;
                        filterLine = (debug) ? filterLine_sse2<uint16_t, 4, 1023, true> : filterLine_sse2<uint16_t, 4, 1023, false>;
                        break;
                    }
                    case 12:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_sse2<uint16_t, true, 4, 4095, true> : filterEdge_sse2<uint16_t, true, 4, 4095, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_sse2<uint16_t, false, 4, 4095, true> : filterEdge_sse2<uint16_t, false, 4, 4095, false>;
                        filterLine = (debug) ? filterLine_sse2<uint16_t, 4, 4095, true> : filterLine_sse2<uint16_t, 4, 4095, false>;
                        break;
                    }
                    case 14:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_sse2<uint16_t, true, 4, 16383, true> : filterEdge_sse2<uint16_t, true, 4, 16383, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_sse2<uint16_t, false, 4, 16383, true> : filterEdge_sse2<uint16_t, false, 4, 16383, false>;
                        filterLine = (debug) ? filterLine_sse2<uint16_t, 4, 16383, true> : filterLine_sse2<uint16_t, 4, 16383, false>;
                        break;
                    }
                    default:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_sse2<uint16_t, true, 4, 65535, true> : filterEdge_sse2<uint16_t, true, 4, 65535, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_sse2<uint16_t, false, 4, 65535, true> : filterEdge_sse2<uint16_t, false, 4, 65535, false>;
                        filterLine = (debug) ? filterLine_sse2<uint16_t, 4, 65535, true> : filterLine_sse2 < uint16_t, 4, 65535, false >;
                        break;
                    }
                }

                filter_ = &BWDIF::filter<uint16_t>;
            }
            break;
            default:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_sse2<float, true, 4, 1, true> : filterEdge_sse2<float, true, 4, 1, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_sse2<float, false, 4, 1, true> : filterEdge_sse2<float, false, 4, 1, false>;
                filterLine = (debug) ? filterLine_sse2<float, 4, 1, true> : filterLine_sse2<float, 4, 1, false>;
                filter_ = &BWDIF::filter<float>;
            }
            break;
        }
    }
    else
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_c<uint8_t, true, 1, 255, true> : filterEdge_c<uint8_t, true, 1, 255, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_c<uint8_t, false, 1, 255, true> : filterEdge_c<uint8_t, false, 1, 255, false>;
                filterLine = (debug) ? filterLine_c<uint8_t, 1, 255, true> : filterLine_c<uint8_t, 1, 255, false>;
                filter_ = &BWDIF::filter<uint8_t>;
            }
            break;
            case 2:
            {
                switch (vi.BitsPerComponent())
                {
                    case 10:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_c<uint16_t, true, 1, 1023, true> : filterEdge_c<uint16_t, true, 1, 1023, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_c<uint16_t, false, 1, 1023, true> : filterEdge_c<uint16_t, false, 1, 1023, false>;
                        filterLine = (debug) ? filterLine_c<uint16_t, 1, 1023, true> : filterLine_c<uint16_t, 1, 1023, false>;
                        break;
                    }
                    case 12:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_c<uint16_t, true, 1, 4095, true> : filterEdge_c<uint16_t, true, 1, 4095, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_c<uint16_t, false, 1, 4095, true> : filterEdge_c<uint16_t, false, 1, 4095, false>;
                        filterLine = (debug) ? filterLine_c<uint16_t, 1, 4095, true> : filterLine_c<uint16_t, 1, 4095, false>;
                        break;
                    }
                    case 14:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_c<uint16_t, true, 1, 16383, true> : filterEdge_c<uint16_t, true, 1, 16383, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_c<uint16_t, false, 1, 16383, true> : filterEdge_c<uint16_t, false, 1, 16383, false>;
                        filterLine = (debug) ? filterLine_c<uint16_t, 1, 16383, true> : filterLine_c<uint16_t, 1, 16383, false>;
                        break;
                    }
                    default:
                    {
                        filterEdgeWithSpat = (debug) ? filterEdge_c<uint16_t, true, 1, 65535, true> : filterEdge_c<uint16_t, true, 1, 65535, false>;
                        filterEdgeWithoutSpat = (debug) ? filterEdge_c<uint16_t, false, 1, 65535, true> : filterEdge_c<uint16_t, false, 1, 65535, false>;
                        filterLine = (debug) ? filterLine_c<uint16_t, 1, 65535, true> : filterLine_c<uint16_t, 1, 65535, false>;
                        break;
                    }
                }

                filter_ = &BWDIF::filter<uint16_t>;
            }
            break;
            default:
            {
                filterEdgeWithSpat = (debug) ? filterEdge_c<float, true, 1, 1, true> : filterEdge_c<float, true, 1, 1, false>;
                filterEdgeWithoutSpat = (debug) ? filterEdge_c<float, false, 1, 1, true> : filterEdge_c<float, false, 1, 1, false>;
                filterLine = (debug) ? filterLine_c<float, 1, 1, true> : filterLine_c<float, 1, 1, false>;
                filter_ = &BWDIF::filter<float>;
            }
            break;
        }
    }

    if (field_ == -2 || field_ > 1)
    {
        vi.num_frames <<= 1;
        vi.fps_numerator <<= 1;
    }

    has_at_least_v8 = env->FunctionExists("propShow");
}

PVideoFrame __stdcall BWDIF::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame edeint{ edeint_ ? edeint_->GetFrame(n, env) : nullptr };

    const int field_no_prop{ [&]()
    {
        if (field_ == -1)
            return child->GetParity(n) ? 1 : 0;
        else if (field_ == -2)
            return child->GetParity(n >> 1) ? 3 : 2;
        else
            return -1;
    } () };

    int field{ (field_ > -1) ? field_ : field_no_prop };
    int tff{ field };

    const int n_orig{ n };
    if (field > 1)
        n >>= 1;

    PVideoFrame prev{ child->GetFrame(std::max(n - 1, 0), env) };
    PVideoFrame cur{ child->GetFrame(n, env) };
    PVideoFrame next{ child->GetFrame(std::min(n + 1, vi.num_frames - 1), env) };
    PVideoFrame dst{ env->NewVideoFrame(vi) };

    if (field_ < 0)
    {
        if (has_at_least_v8)
        {
            int err;
            const int64_t field_based{ env->propGetInt(env->getFramePropsRO(cur), "_FieldBased", 0, &err) };
            if (err == 0)
            {
                switch (field_based)
                {
                    case 1:
                    {
                        field = 0;
                        tff = 0;
                        break;
                    }
                    case 2:
                    {
                        field = 1;
                        tff = 1;
                        break;
                    }
                    default:
                    {
                        if (pass_)
                            return cur;
                        else
                            env->ThrowError("BWDIF: _FieldBased frame property must be greater than 0.");
                        break;
                    }
                }

                if (field_ > 1 || field_no_prop > 1)
                {
                    if (field_based == 0)
                        field -= 2;

                    field = (n_orig & 1) ? (1 - field) : field;
                }
            }
            else
            {
                if (field > 1)
                {
                    field -= 2;
                    field = (n_orig & 1) ? (1 - field) : field;
                }
            }
        }
        else
        {
            if (field > 1)
            {
                field -= 2;
                field = (n_orig & 1) ? (1 - field) : field;
            }
        }
    }
    else
    {
        if (field > 1)
        {
            field -= 2;
            field = (n_orig & 1) ? (1 - field) : field;
        }
    }

    if (tff > 1)
        tff = (tff == 2) ? 1 : 0;
    else
        tff = 1 - tff;

    (this->*filter_)(prev, cur, next, dst, edeint, field, tff, env);

    if (has_at_least_v8)
    {
        env->copyFrameProps(cur, dst);
        env->propSetInt(env->getFramePropsRW(dst), "_FieldBased", 0, 0);

        if (field_ > 1 || field_no_prop > 1)
        {
            AVSMap* props{ env->getFramePropsRW(dst) };
            int errNum;
            int errDen;
            int64_t durationNum{ env->propGetInt(props, "_DurationNum", 0, &errNum) };
            int64_t durationDen{ env->propGetInt(props, "_DurationDen", 0, &errDen) };
            if (errNum == 0 && errDen == 0)
            {
                muldivRational(&durationNum, &durationDen, 1, 2);
                env->propSetInt(props, "_DurationNum", durationNum, 0);
                env->propSetInt(props, "_DurationNum", durationNum, 0);
            }
        }
    }

    return dst;
}

AVSValue __cdecl Create_BWDIF(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    enum { Clip, Field, Edeint, Opt, Thr, Debug, Pass };

    auto edeint{ (args[Edeint].Defined()) ? args[Edeint].AsClip() : nullptr };

    return new BWDIF(
        args[Clip].AsClip(),
        args[Field].AsInt(-1),
        edeint,
        args[Opt].AsInt(-1),
        args[Thr].AsFloatf(0.0f),
        args[Debug].AsBool(false),
        args[Pass].AsBool(false),
        env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("BWDIF", "c[field]i[edeint]c[opt]i[thr]f[debug]b[pass]b", Create_BWDIF, 0);

    return "BWDIF";
}
