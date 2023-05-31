#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "Arduino.h"

const float MIN_PRECISION = 2.0 * __FLT_EPSILON__;

// simpson's rule
const double a2_0 = 1.0 / 6.0;
const double a2_1 = 2.0 / 3.0;
// const double a2_2 = a2_0;

// 3-step newton-cotes
const double a3_0 = 1.0 / 8.0;
const double a3_1 = 3.0 / 8.0;
// const double a3_2 = a3_1;
// const double a3_3 = a3_0;

// 4-step method
const double a4_0 = 11.0 / 120.0;
const double a4_1 = 27.0 / 40.0;
const double a4_2 = -8.0 / 15.0;
// const double a4_3 = a4_1;
// const double a4_4 = a4_0;

// 5-step method
const double a5_0 = 151.0 / 4320.0;
const double a5_1 = 539.0 / 2160.0;
const double a5_2 = 931.0 / 4320.0;
// const double a5_3 = a5_2;
// const double a5_4 = a5_1;
// const double a5_5 = a5_0;

// 6-step method
const double a6_0 = 31.0 / 720.0;
const double a6_1 = 16807.0 / 79200.0;
const double a6_2 = 243.0 / 1760.0;
const double a6_3 = 16.0 / 75.0;
// const double a6_4 = a6_2;
// const double a6_5 = a6_1;
// const double a6_6 = a6_0;

// 7-step method
const double a7_0 = 739.0 / 17280.0;
const double a7_1 = 40817.0 / 190080.0;
const double a7_2 = 729.0 / 7040.0;
const double a7_3 = 2401.0 / 17280.0;
// const double a7_4 = a7_3;
// const double a7_5 = a7_2;
// const double a7_6 = a7_1;
// const double a7_7 = a7_0;


// 4th-order adaptive quadrature
template <typename Function>
float adapt4(Function&& f, float a, float b, float precision) {
    float result = 0;
    float h = b - a;
    float x_0 = a;
    double f_start = f(a);


    while (x_0 != b) {
        // step to next interval
        float x_1 = x_0 + h;
        if (x_1 > b) {
            x_1 = b;
        }
        float error_target = h * precision / (b - a);
        // evaluate function within interval
        double f_1_3 = f(x_0 + h / 3.0);
        double f_1_2 = f(x_0 + h / 2.0);
        double f_2_3 = f(x_0 + 2.0 * h / 3.0);
        double f_end = f(x_1);
        // make 3-point and 2-point estimates
        double s3 = h * (a3_0 * f_start + a3_1 * f_1_3 + a3_1 * f_2_3 + a3_0 * f_end);
        double s2 = h * (a2_0 * f_start + a2_1 * f_1_2 + a2_0 * f_end);
        // update interval size based on estimated error
        float error = 0.8 * abs(s3 - s2);
        if (error == 0) {
            h = b - x_0;
        }
        else {
            h *= 0.9 * sqrt(sqrt(error_target / error));
        }
        if (error < error_target) {
            result += s3;
            x_0 = x_1;
            f_start = f_end;
        }
    }

    return result;
}

template <typename Function>
float vsvoq2(Function&& f, float a, float b, float precision) {
    float result = 0;
    float h = b - a;
    float x_0 = a;
    double f_start = f(a);

    uint8_t next_order = 4;

    // keep track of how many iterations we stay before moving on
    // if you've been there too long, prohibit higher-order approaches
    size_t steps_here = 0;

    while (x_0 != b) {
        float x_1 = x_0 + h;
        if (x_1 > b) {
            h = b - x_0;
            x_1 = b;
        }
        float error_target = h * precision / (b - a);
        float error;
        float approx;

        double f_1_7 = f(x_0 + h / 7.0);
        double f_1_3 = f(x_0 + h / 3.0);
        double f_3_7 = f(x_0 + 3.0 * h / 7.0);
        double f_1_2 = f(x_0 + h / 2.0);
        double f_4_7 = f(x_0 + 4.0 * h / 7.0);
        double f_2_3 = f(x_0 + 2.0 * h / 3.0);
        double f_6_7 = f(x_0 + 6.0 * h / 7.0);
        double f_end = f(x_1);

        double h_8;
        double h_6;
        double h_4;
        double h_2;

        // compute estimates for each of 4 approximations
        if (steps_here <= 1) {
            double s7 = h * (
                a7_0 * f_start + a7_1 * f_1_7 + a7_2 * f_1_3 + a7_3 * f_3_7
                + a7_3 * f_4_7 + a7_2 * f_2_3 + a7_1 * f_6_7 + a7_0 * f_end
            );
            double s6 = h * (
                a6_0 * f_start + a6_1 * f_1_7 + a6_2 * f_1_3 + a6_3 * f_1_2
                + a6_2 * f_2_3 + a6_1 * f_6_7 + a6_0 * f_end
            );

            // estimated error for 8th-order method
            float error_8 = abs((72.0 / 5.0) * (s7 - s6));
            // next stepsize for 8th-order method
            h_8 = (error_8 == 0) ? b - x_0 : h * 0.9 * pow(error_target / error_8, 1.0 / 8.0);

            if (next_order == 8) {
                error = error_8;
                approx = s7;
            }
        }

        if (steps_here <= 3) {
            double s5 = h * (
                a5_0 * f_start + a5_1 * f_1_7 + a5_2 * f_3_7
                + a5_2 * f_4_7 + a5_1 * f_6_7 + a5_0 * f_end
            );
            double s4 = h * (
                a4_0 * f_start + a4_1 * f_1_3 + a4_2 * f_1_2
                + a4_1 * f_2_3 + a4_0 * f_end
            );

            // estimated error for 6th-order method
            float error_6 = abs((54.0 / 397.0) * (s5 - s4));
            // next stepsize for 6th-order method
            h_6 = (error_6 == 0) ? b - x_0 : h * 0.9 * pow(error_target / error_6, 1.0 / 6.0);

            if (next_order == 6) {
                error = error_6;
                approx = s5;
            }

            double s3 = h * (
                a3_0 * f_start + a3_1 * f_1_3 + a3_1 * f_2_3 + a3_0 * f_end
            );
            double s2 = h * (
                a2_0 * f_start + a2_1 * f_1_2 + a2_0 * f_end
            );

            // estimated error for 4th-order method
            float error_4 = abs((4.0 / 5.0) * (s3 - s2));
            // next stepsize for 4th-order method
            h_4 = (error_4 == 0) ? b - x_0 : h * 0.9 * sqrt(sqrt(error_target / error_4));

            if (next_order == 4) {
                error = error_4;
                approx = s3;
            }
        }
        
        double s1 = h * (0.5 * f_start + 0.5 * f_end);
        double s0 = h * f_1_2;

        // estimated error for 2nd-order method
        float error_2 = abs(s1 - s0) / 3.0;
        // next stepsize for 2nd-order method
        h_2 = (error_2 == 0) ? b - x_0 : h * 0.9 * sqrt(error_target / error_2);

        if (next_order == 2) {
            error = error_2;
            approx = s1;
        }

        // if error for current order is within tolerance, use it
        if (error <= error_target) {
            result += approx;
            x_0 = x_1;
            f_start = f_end;
            steps_here = 0;
        }
        else {
            steps_here++;
        }

        // select next order to maximize step size
        next_order = 2;
        h = h_2;
        if (steps_here > 3) {
            continue;
        }
        if (h_4 > h) {
            next_order = 4;
            h = h_4;
        }
        if (h_6 > h) {
            next_order = 6;
            h = h_6;
        }
        if (steps_here > 1) {
            continue;
        }
        if (h_8 > h_2) {
            next_order = 8;
            h = h_8;
        }
    }

    return result;
}

template <typename Function>
float vsvoq3(Function&& f, float a, float b, float precision) {
    return 0;
}

template <typename Function>
float integrate(Function&& f, float a, float b, float precision) {
    precision = max(precision, MIN_PRECISION);
    // decide whether to hand off to VSVO, or just adaptive 4th-order rule
    if (precision >= 1e-4) {
        return adapt4(f, a, b, precision);
    }
    else {
        return vsvoq2(f, a, b, precision);
    }
}

#endif