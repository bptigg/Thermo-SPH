#pragma once
#include <cuda_runtime.h>
#include "kernel.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// Inline GPU device kernel evaluator (Zero vtable overhead)
struct GPUKernelEvaluator {
    KernelType type;

    __device__ inline float value(float dist, float h) const {
        float q = dist / h;
        if (q >= 2.0f) return 0.0f;

        if (type == KernelType::CubicSpline) {
            float sigma = 10.0f / (7.0f * static_cast<float>(M_PI) * h * h);
            if (q < 1.0f) {
                return sigma * (1.0f - 1.5f * q * q + 0.75f * q * q * q);
            }
            float term = 2.0f - q;
            return sigma * 0.25f * term * term * term;
        } 
        else { // WendlandC2
            float sigma = 7.0f / (4.0f * static_cast<float>(M_PI) * h * h);
            float term = 1.0f - 0.5f * q;
            return sigma * powf(term, 4.0f) * (1.0f + 2.0f * q);
        }
    }

    __device__ inline void gradient(float rx, float ry, float dist, float h, float& gradX, float& gradY) const {
        gradX = 0.0f;
        gradY = 0.0f;
        if (dist <= 1e-6f || dist >= 2.0f * h) return;

        float q = dist / h;
        float dWdq = 0.0f;

        if (type == KernelType::CubicSpline) {
            float sigma = 10.0f / (7.0f * static_cast<float>(M_PI) * h * h);
            dWdq = (q < 1.0f) ? sigma * (-3.0f * q + 2.25f * q * q)
                              : sigma * (-0.75f * (2.0f - q) * (2.0f - q));
        } 
        else { // WendlandC2
            float sigma = 7.0f / (4.0f * static_cast<float>(M_PI) * h * h);
            float term = 1.0f - 0.5f * q;
            dWdq = sigma * (-5.0f * q * powf(term, 3.0f));
        }

        float factor = dWdq / (h * dist);
        gradX = rx * factor;
        gradY = ry * factor;
    }
};