#include <iostream>
#include <vector>
#include <format>

typedef double (*Function)(double, double);

double Func(double x, double y)
{
	return y - sin(x * x * x * x);
}

bool Refine(double* x0, double* y0, Function func)
{
	static constexpr double step = 1e-9;

	static std::vector<double> zs{};
	static std::vector<double> dels{};
	static std::vector<double> qs{};
	zs.clear(); dels.clear(); qs.clear();

	double& xk = *x0, & yk = *y0;
	int k;
	for (k = 0;; k++)
	{
		// Evaluate function
		double z = func(xk, yk);
		if (!std::isfinite(z)) return false;

		zs.push_back(z);
		if (z == 0) break;

		double zxn = func(xk + xk * step, yk);
		double zyn = func(xk, yk + yk * step);

		if (zxn == z && zyn == z) break; // Gradient is zero, newton's method will fail
		if (!std::isfinite(zxn) || !std::isfinite(zyn)) break;

		// Estimate partial derivatives
		double dx = (zxn - z) / (xk * step);
		double dy = (zyn - z) / (yk * step);

		// Compute values
		double absD = dx * dx + dy * dy;
		double fac = z / absD;

		double delX = dx * fac;
		double delY = dy * fac;

		// Apply deltas
		xk -= delX;
		yk -= delY;

		// Calculate distance
		double xErr = abs((nextafter(xk, DBL_MAX) - xk) / delX);
		double yErr = abs((nextafter(yk, DBL_MAX) - yk) / delY);
		if (xErr > 0.01 && yErr > 0.01) break; // Step is below relative threshold
		dels.push_back(abs(z) / sqrt(absD));

		// Add q
		if (k >= 2)
			qs.push_back(log(abs(dels[k] / dels[k - 1])) / log(abs(dels[k - 1] / dels[k - 2])));

		if (k == 1199) break;
	}
	
	// Process information to determine convergence
	double descentRate = std::max(log(zs[k - 1]) - log(zs[k]), log(zs[k - 2]) - log(zs[k - 1]));
}

int main()
{
	double x = 12.1, y = 10.1;

	if (Refine(&x, &y, Func))
	{
		std::cout << std::format("Finished: ({:#.17g}, {:#.17g})\n", x, y);
	}
	else
	{
		std::cout << std::format("Failed to converge: ({:#.17g}, {:#.17g})\n", x, y);
	}
}