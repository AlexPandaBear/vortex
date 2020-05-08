#include <iostream>
#include "SimManager.hxx"

int main(int argc, char const *argv[])
{
	SimManager SM;

	for (size_t i = 0; i < 10; i++)
	{
		for (size_t j = -5; j <= 5; j++)
		{
			SM.addVtx(i, j, exp(-j*j), 0.01);
		}
	}

	SM.buildTimeSample(0., 10., 20);
	SM.setXPeriodicityTo(true, 10.);
	SM.setMethodTo("rk4");
	SM.sim();

	std::cout << SM.computeVorticityAt(2, 5, 0, 1., true, 10.) << std::endl;

	return 0;
}