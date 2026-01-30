#include "roots.hpp"
#include <cmath>

constexpr int MAX_ITER = 100;
constexpr double TOL = 1e-6;

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double mid = (a + b) / 2.0;
        double fmid = f(mid);

        if (std::fabs(fmid) < TOL)
        {
            *root = mid;
            return true;
        }

        if (fa * fmid < 0)
        {
            b = mid;
            fb = fmid;
        }
        else
        {
            a = mid;
            fa = fmid;
        }
    }

    *root = (a + b) / 2.0;
    return true;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double c = b - fb * (b - a) / (fb - fa);
        double fc = f(c);

        if (std::fabs(fc) < TOL)
        {
            *root = c;
            return true;
        }

        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
            fa *= 0.5;   
        }
        else
        {
            a = c;
            fa = fc;
            fb *= 0.5;   
        }
    }

    *root = b - fb * (b - a) / (fb - fa);
    return true;
}


bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    double x = c;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double fx = f(x);
        double gx = g(x);

        if (std::fabs(gx) < TOL)
            return false;

        double x_next = x - fx / gx;

        if (x_next < a || x_next > b)
            return false;

        if (std::fabs(f(x_next)) < TOL)
        {
            *root = x_next;
            return true;
        }

        x = x_next;
    }

    *root = x;
    return true;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    double x0 = c;
    double x1 = b;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < MAX_ITER; ++i)
    {
        if (std::fabs(f1 - f0) < TOL)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        if (std::fabs(f(x2)) < TOL)
        {
            *root = x2;
            return true;
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    *root = x1;
    return true;
}
