#include <iostream>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char *argv[])
{

  const double pressure = 1 * AtmosphericPressure;
  const double temperature = 293.15;

  // Setup the gas.
  MediumMagboltz gas("Ar", 89., "CO2", 10., "CH4", 1.);
  gas.SetTemperature(temperature);
  gas.SetPressure(pressure);

  // Set the field range to be covered by the gas table.
  const size_t nE = 40;
  const double emin = 100.;
  const double emax = 100000.;
  // Flag to request logarithmic spacing.
  constexpr bool useLog = true;
  gas.SetFieldGrid(emin, emax, nE, useLog, 0, 3, 6);

  const int ncoll = 10;
  // Run Magboltz to generate the gas table.
  gas.GenerateGasTable(ncoll);
  // Save the table.
  gas.WriteGasFile("Gasfile/ar_89_co2_10_ch4_1_3T.gas");
}
