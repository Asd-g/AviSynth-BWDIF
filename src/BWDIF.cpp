#include <algorithm>

#include "avisynth.h"
#include "BWDIF.h"

template<typename pixel_t, bool spat>
static inline void filterEdge_c(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint_, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept
{
    const pixel_t* prev2 = reinterpret_cast<const pixel_t*>(_prev2);
    const pixel_t* prev = reinterpret_cast<const pixel_t*>(_prev);
    const pixel_t* cur = reinterpret_cast<const pixel_t*>(_cur);
    const pixel_t* edeint = reinterpret_cast<const pixel_t*>(edeint_);
    const pixel_t* next = reinterpret_cast<const pixel_t*>(_next);
    const pixel_t* next2 = reinterpret_cast<const pixel_t*>(_next2);
    pixel_t* __restrict dst = reinterpret_cast<pixel_t*>(_dst);

    const pixel_t* prev2Above2 = prev2 - stride2;
    const pixel_t* prev2Below2 = prev2 + stride2;
    const pixel_t* prevAbove = prev + negativeStride;
    const pixel_t* prevBelow = prev + positiveStride;
    const pixel_t* curAbove = cur + negativeStride;
    const pixel_t* curBelow = cur + positiveStride;
    const pixel_t* nextAbove = next + negativeStride;
    const pixel_t* nextBelow = next + positiveStride;
    const pixel_t* next2Above2 = next2 - stride2;
    const pixel_t* next2Below2 = next2 + stride2;

    for (int x = 0; x < width; x++)
    {
        if constexpr (std::is_integral_v<pixel_t>)
        {
            const int c = curAbove[x];
            const int d = (prev2[x] + next2[x]) >> 1;
            const int e = curBelow[x];
            const int temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const int temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) >> 1;
            const int temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) >> 1;
            int diff = std::max({ temporal_diff0 >> 1, temporal_diff1, temporal_diff2 });

            if (!diff)
                dst[x] = d;
            else
            {
                if constexpr (spat)
                {
                    const int b = ((prev2Above2[x] + next2Above2[x]) >> 1) - c;
                    const int f = ((prev2Below2[x] + next2Below2[x]) >> 1) - e;
                    const int dc = d - c;
                    const int de = d - e;
                    const int maximum = std::max({ de, dc, std::min(b, f) });
                    const int minimum = std::min({ de, dc, std::max(b, f) });
                    diff = std::max({ diff, minimum, -maximum });
                }

                int interpol = std::clamp(!edeint ? (c + e) >> 1 : static_cast<int>(edeint[x]), d - diff, d + diff);
                dst[x] = std::clamp(interpol, 0, peak);
            }
        }
        else
        {
            const float c = curAbove[x];
            const float d = (prev2[x] + next2[x]) * 0.5f;
            const float e = curBelow[x];
            const float temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const float temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) * 0.5f;
            const float temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) * 0.5f;
            float diff = std::max({ temporal_diff0 * 0.5f, temporal_diff1, temporal_diff2 });

            if (!diff)
                dst[x] = d;
            else
            {
                if constexpr (spat)
                {
                    const float b = ((prev2Above2[x] + next2Above2[x]) * 0.5f) - c;
                    const float f = ((prev2Below2[x] + next2Below2[x]) * 0.5f) - e;
                    const float dc = d - c;
                    const float de = d - e;
                    const float maximum = std::max({ de, dc, std::min(b, f) });
                    const float minimum = std::min({ de, dc, std::max(b, f) });
                    diff = std::max({ diff, minimum, -maximum });
                }

                dst[x] = std::clamp(!edeint ? (c + e) * 0.5f : edeint[x], d - diff, d + diff);
            }
        }
    }
}

template<typename pixel_t>
static inline void filterLine_c(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint_, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept
{
    const pixel_t* prev2 = reinterpret_cast<const pixel_t*>(_prev2);
    const pixel_t* prev = reinterpret_cast<const pixel_t*>(_prev);
    const pixel_t* cur = reinterpret_cast<const pixel_t*>(_cur);
    const pixel_t* edeint = reinterpret_cast<const pixel_t*>(edeint_);
    const pixel_t* next = reinterpret_cast<const pixel_t*>(_next);
    const pixel_t* next2 = reinterpret_cast<const pixel_t*>(_next2);
    
    pixel_t* __restrict dst = reinterpret_cast<pixel_t*>(_dst);

    const pixel_t* prev2Above4 = prev2 - stride4;
    const pixel_t* prev2Above2 = prev2 - stride2;
    const pixel_t* prev2Below2 = prev2 + stride2;
    const pixel_t* prev2Below4 = prev2 + stride4;
    const pixel_t* prevAbove = prev - stride;
    const pixel_t* prevBelow = prev + stride;
    const pixel_t* curAbove3 = cur - stride3;
    const pixel_t* curAbove = cur - stride;
    const pixel_t* curBelow = cur + stride;
    const pixel_t* curBelow3 = cur + stride3;
    const pixel_t* nextAbove = next - stride;
    const pixel_t* nextBelow = next + stride;
    const pixel_t* next2Above4 = next2 - stride4;
    const pixel_t* next2Above2 = next2 - stride2;
    const pixel_t* next2Below2 = next2 + stride2;
    const pixel_t* next2Below4 = next2 + stride4;

    for (int x = 0; x < width; x++)
    {
        if constexpr (std::is_integral_v<pixel_t>)
        {
            const int c = curAbove[x];
            const int d = (prev2[x] + next2[x]) >> 1;
            const int e = curBelow[x];
            const int temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const int temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) >> 1;
            const int temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) >> 1;
            int diff = std::max({ temporal_diff0 >> 1, temporal_diff1, temporal_diff2 });

            if (!diff)
                dst[x] = d;
            else
            {
                const int b = ((prev2Above2[x] + next2Above2[x]) >> 1) - c;
                const int f = ((prev2Below2[x] + next2Below2[x]) >> 1) - e;
                const int dc = d - c;
                const int de = d - e;
                const int maximum = std::max({ de, dc, std::min(b, f) });
                const int minimum = std::min({ de, dc, std::max(b, f) });
                diff = std::max({ diff, minimum, -maximum });

                int interpol;
                if (!edeint)
                {
                    if (std::abs(c - e) > temporal_diff0)
                        interpol = (((coef_hf[0] * (prev2[x] + next2[x])
                            - coef_hf[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                            + coef_hf[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) >> 2)
                            + coef_lf[0] * (c + e) - coef_lf[1] * (curAbove3[x] + curBelow3[x])) >> 13;
                    else
                        interpol = (coef_sp[0] * (c + e) - coef_sp[1] * (curAbove3[x] + curBelow3[x])) >> 13;

                    interpol = std::clamp(interpol, d - diff, d + diff);
                }
                else
                    interpol = std::clamp(static_cast<int>(edeint[x]), d - diff, d + diff);

                dst[x] = std::clamp(interpol, 0, peak);
            }
        }
        else {
            const float c = curAbove[x];
            const float d = (prev2[x] + next2[x]) * 0.5f;
            const float e = curBelow[x];
            const float temporal_diff0 = std::abs(prev2[x] - next2[x]);
            const float temporal_diff1 = (std::abs(prevAbove[x] - c) + std::abs(prevBelow[x] - e)) * 0.5f;
            const float temporal_diff2 = (std::abs(nextAbove[x] - c) + std::abs(nextBelow[x] - e)) * 0.5f;
            float diff = std::max({ temporal_diff0 * 0.5f, temporal_diff1, temporal_diff2 });

            if (!diff)
                dst[x] = d;
            else
            {
                const float b = ((prev2Above2[x] + next2Above2[x]) * 0.5f) - c;
                const float f = ((prev2Below2[x] + next2Below2[x]) * 0.5f) - e;
                const float dc = d - c;
                const float de = d - e;
                const float maximum = std::max({ de, dc, std::min(b, f) });
                const float minimum = std::min({ de, dc, std::max(b, f) });
                diff = std::max({ diff, minimum, -maximum });

                float interpol;
                if (std::abs(c - e) > temporal_diff0)
                    interpol = ((coef_hf_f[0] * (prev2[x] + next2[x])
                        - coef_hf_f[1] * (prev2Above2[x] + next2Above2[x] + prev2Below2[x] + next2Below2[x])
                        + coef_hf_f[2] * (prev2Above4[x] + next2Above4[x] + prev2Below4[x] + next2Below4[x])) * 0.25f
                        + coef_lf_f[0] * (c + e) - coef_lf_f[1] * (curAbove3[x] + curBelow3[x]));
                else
                    interpol = coef_sp_f[0] * (c + e) - coef_sp_f[1] * (curAbove3[x] + curBelow3[x]);

                dst[x] = std::clamp(!edeint ? interpol : edeint[x], d - diff, d + diff);
            }
        }
    }
}

class BWDIF : public GenericVideoFilter
{
    int field_;
    PClip edeint_;
    int opt_;
    int edgeStep, lineStep, peak;
    bool has_at_least_v8;

    void (*filterEdgeWithSpat)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
    void (*filterEdgeWithoutSpat)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int positiveStride, const int negativeStride, const int stride2, const int step, const int peak) noexcept;
    void (*filterLine)(const void* _prev2, const void* _prev, const void* _cur, const void* _next, const void* _next2, void* _dst, const void* edeint, const int width, const int stride, const int stride2, const int stride3, const int stride4, const int step, const int peak) noexcept;
    template<typename pixel_t>
    void filter(PVideoFrame& prevFrame, PVideoFrame& curFrame, PVideoFrame& nextFrame, PVideoFrame& dstFrame, PVideoFrame& edeint, const int field, const BWDIF* const __restrict, IScriptEnvironment* env) noexcept;

public:
    BWDIF(PClip _child, int field, PClip edeint, int opt, IScriptEnvironment* env);
    PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
    int __stdcall SetCacheHints(int cachehints, int frame_range)
    {
        return cachehints == CACHE_GET_MTMODE ? MT_NICE_FILTER : 0;
    }
};

template<typename pixel_t>
void BWDIF::filter(PVideoFrame& prevFrame, PVideoFrame& curFrame, PVideoFrame& nextFrame, PVideoFrame& dstFrame, PVideoFrame& edeintFrame, const int field, const BWDIF* const __restrict, IScriptEnvironment* env) noexcept
{
    int planes_y[4] = { PLANAR_Y, PLANAR_U, PLANAR_V, PLANAR_A };
    int planes_r[4] = { PLANAR_G, PLANAR_B, PLANAR_R, PLANAR_A };
    const int* current_planes = (vi.IsYUV() || vi.IsYUVA()) ? planes_y : planes_r;
    const int planecount = std::min(vi.NumComponents(), 3);
    for (int i = 0; i < planecount; i++)
    {
        const int plane = current_planes[i];

        const int stride = curFrame->GetPitch(plane) / sizeof(pixel_t);
        const int edeint_stride = edeintFrame ? edeintFrame->GetPitch(plane) / sizeof(pixel_t) : 0;
        const int dst_stride = dstFrame->GetPitch(plane) / sizeof(pixel_t);
        const int width = curFrame->GetRowSize(plane) / sizeof(pixel_t);
        const int height = curFrame->GetHeight(plane);
        const pixel_t* prev = reinterpret_cast<const pixel_t*>(prevFrame->GetReadPtr(plane));
        const pixel_t* cur = reinterpret_cast<const pixel_t*>(curFrame->GetReadPtr(plane));
        const pixel_t* next = reinterpret_cast<const pixel_t*>(nextFrame->GetReadPtr(plane));
        const pixel_t* edeint = edeintFrame ? reinterpret_cast<const pixel_t*>(edeintFrame->GetReadPtr(plane)) : 0;
        pixel_t* __restrict dst = reinterpret_cast<pixel_t*>(dstFrame->GetWritePtr(plane));

        env->BitBlt(dstFrame->GetWritePtr(plane) + dstFrame->GetPitch(plane) * (static_cast<int64_t>(1) - field), dstFrame->GetPitch(plane) * 2, curFrame->GetReadPtr(plane) + curFrame->GetPitch(plane) * (static_cast<int64_t>(1) - field), curFrame->GetPitch(plane) * 2,
            curFrame->GetRowSize(plane), height / 2);

        prev += static_cast<int64_t>(stride) * field;
        cur += static_cast<int64_t>(stride) * field;
        next += static_cast<int64_t>(stride) * field;
        edeint += static_cast<int64_t>(edeint_stride) * field;
        dst += static_cast<int64_t>(dst_stride) * field;

        const pixel_t* prev2 = field ? prev : cur;
        const pixel_t* next2 = field ? cur : next;

        for (int y = field; y < height; y += 2) {
            if ((y < 4) || (y + 5 > height)) {
                if ((y < 2) || (y + 3 > height))
                    filterEdgeWithoutSpat(prev2, prev, cur, next, next2, dst, edeint, width,
                        y + 1 < height ? stride : -stride,
                        y > 0 ? -stride : stride,
                        stride * 2,
                        edgeStep, peak);
                else
                    filterEdgeWithSpat(prev2, prev, cur, next, next2, dst, edeint, width,
                        y + 1 < height ? stride : -stride,
                        y > 0 ? -stride : stride,
                        stride * 2,
                        edgeStep, peak);
            }
            else {
                filterLine(prev2, prev, cur, next, next2, dst, edeint, width,
                    stride, stride * 2, stride * 3, stride * 4,
                    lineStep, peak);
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

BWDIF::BWDIF(PClip _child, int field, PClip edeint, int opt, IScriptEnvironment* env)
    : GenericVideoFilter(_child), field_(field), edeint_(edeint), opt_(opt)
{
    if (!vi.IsPlanar())
        env->ThrowError("BWDIF: only planar formats are supported.");
    if (vi.height < 4)
        env->ThrowError("BWDIF: height must be greater than or equal to 4.");
    if (field_ < -2 || field_ > 3)
        env->ThrowError("BWDIF: field must be -2, -1, 0, 1, 2, or 3.");
    if (opt_ < -1 || opt_ > 3)
        env->ThrowError("BWDIF: opt must be between -1..3.");
    if (!(env->GetCPUFlags() & CPUF_AVX512F) && opt_ == 3)
        env->ThrowError("BWDIF: opt=3 requires AVX512F.");
    if (!(env->GetCPUFlags() & CPUF_AVX2) && opt_ == 2)
        env->ThrowError("BWDIF: opt=2 requires AVX2.");
    if (!(env->GetCPUFlags() & CPUF_SSE2) && opt_ == 1)
        env->ThrowError("BWDIF: opt=1 requires SSE2.");

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

    edgeStep = lineStep = 0;

    if ((!!(env->GetCPUFlags() & CPUF_AVX512F) && opt_ < 0) || opt_ == 3)
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = filterEdge_avx512<uint8_t, true>;
                filterEdgeWithoutSpat = filterEdge_avx512<uint8_t, false>;
                filterLine = filterLine_avx512<uint8_t>;
                edgeStep = 32;
            }
            break;
            case 2:
            {
                filterEdgeWithSpat = filterEdge_avx512<uint16_t, true>;
                filterEdgeWithoutSpat = filterEdge_avx512<uint16_t, false>;
                filterLine = filterLine_avx512<uint16_t>;
                edgeStep = 16;
            }
            break;
            default:
            {
                filterEdgeWithSpat = filterEdge_avx512<float, true>;
                filterEdgeWithoutSpat = filterEdge_avx512<float, false>;
                filterLine = filterLine_avx512<float>;
                edgeStep = 16;
            }
            break;
        }

        lineStep = 16;
    }
    else if ((!!(env->GetCPUFlags() & CPUF_AVX2) && opt_ < 0) || opt_ == 2)
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = filterEdge_avx2<uint8_t, true>;
                filterEdgeWithoutSpat = filterEdge_avx2<uint8_t, false>;
                filterLine = filterLine_avx2<uint8_t>;
                edgeStep = 16;
            }
            break;
            case 2:
            {
                filterEdgeWithSpat = filterEdge_avx2<uint16_t, true>;
                filterEdgeWithoutSpat = filterEdge_avx2<uint16_t, false>;
                filterLine = filterLine_avx2<uint16_t>;
                edgeStep = 8;
            }
            break;
            default:
            {
                filterEdgeWithSpat = filterEdge_avx2<float, true>;
                filterEdgeWithoutSpat = filterEdge_avx2<float, false>;
                filterLine = filterLine_avx2<float>;
                edgeStep = 8;
            }
            break;
        }

        lineStep = 8;
    }
    else if ((!!(env->GetCPUFlags() & CPUF_SSE2) && opt_ < 0) || opt_ == 1)
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = filterEdge_sse2<uint8_t, true>;
                filterEdgeWithoutSpat = filterEdge_sse2<uint8_t, false>;
                filterLine = filterLine_sse2<uint8_t>;
                edgeStep = 8;
            }
            break;
            case 2:
            {
                filterEdgeWithSpat = filterEdge_sse2<uint16_t, true>;
                filterEdgeWithoutSpat = filterEdge_sse2<uint16_t, false>;
                filterLine = filterLine_sse2<uint16_t>;
                edgeStep = 4;
            }
            break;
            default:
            {
                filterEdgeWithSpat = filterEdge_sse2<float, true>;
                filterEdgeWithoutSpat = filterEdge_sse2<float, false>;
                filterLine = filterLine_sse2<float>;
                edgeStep = 4;
            }
            break;
        }

        lineStep = 4;
    }
    else
    {
        switch (vi.ComponentSize())
        {
            case 1:
            {
                filterEdgeWithSpat = filterEdge_c<uint8_t, true>;
                filterEdgeWithoutSpat = filterEdge_c<uint8_t, false>;
                filterLine = filterLine_c<uint8_t>;
            }
            break;
            case 2:
            {
                filterEdgeWithSpat = filterEdge_c<uint16_t, true>;
                filterEdgeWithoutSpat = filterEdge_c<uint16_t, false>;
                filterLine = filterLine_c<uint16_t>;
            }
            break;
            default:
            {
                filterEdgeWithSpat = filterEdge_c<float, true>;
                filterEdgeWithoutSpat = filterEdge_c<float, false>;
                filterLine = filterLine_c<float>;
            }
            break;
        }
    }

    if (field_ == - 2 || field_ > 1)
    {
        vi.num_frames *= 2;
        vi.fps_numerator *= 2;
    }
    if (field_ == -1)
        field_ = child->GetParity(0) ? 1 : 0;
    if (field_ == -2)
        field_ = child->GetParity(0) ? 3 : 2;

    peak = (1 << vi.BitsPerComponent()) - 1;

    has_at_least_v8 = true;
    try { env->CheckVersion(8); }
    catch (const AvisynthError&) { has_at_least_v8 = false; }
}

PVideoFrame __stdcall BWDIF::GetFrame(int n, IScriptEnvironment* env)
{
    PVideoFrame edeint = edeint_ ? edeint_->GetFrame(n, env) : nullptr;

    int field = field_;
    if (field_ == -2 || field_ > 1)
    {
        field -= 2;
        field = n & 1 ? (field == 0) : (field == 1);
        n /= 2;
    }

    PVideoFrame prev = child->GetFrame(std::max(n - 1, 0), env);
    PVideoFrame cur = child->GetFrame(n, env);
    PVideoFrame next = child->GetFrame(std::min(n + 1, vi.num_frames - 1), env);
    PVideoFrame dst = has_at_least_v8 ? env->NewVideoFrameP(vi, &cur) : env->NewVideoFrame(vi);    

    switch (vi.ComponentSize())
    {
        case 1: filter<uint8_t>(prev, cur, next, dst, edeint, field, 0, env); break;
        case 2: filter<uint16_t>(prev, cur, next, dst, edeint, field, 0, env); break;
        default: filter<float>(prev, cur, next, dst, edeint, field, 0, env); break;
    }

    return dst;
}

AVSValue __cdecl Create_BWDIF(AVSValue args, void* user_data, IScriptEnvironment* env)
{
    return new BWDIF(
        args[0].AsClip(),
        args[1].AsInt(-1),
        args[2].AsClip(),
        args[3].AsInt(-1),
        env);
}

const AVS_Linkage* AVS_linkage;

extern "C" __declspec(dllexport)
const char* __stdcall AvisynthPluginInit3(IScriptEnvironment * env, const AVS_Linkage* const vectors)
{
    AVS_linkage = vectors;

    env->AddFunction("BWDIF", "c[field]i[edeint]c[opt]i", Create_BWDIF, 0);

    return "BWDIF";
}
