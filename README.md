# From-Interdependent-Networks-to-Two-Interaction-Physical-Systems
Simulation codes of the network system employed to obtain the theoretical results reported in the article.


main scripts :

one_network_SC_temperature.m :   Evaluation of the effective resistance in a superconducting square-grid network under constant applied current, as a function of varying heat-bath temperature.

one_network_SC_current.m :  Evaluation of the effective resistance in a superconducting square-grid network under constant applied heat-bath temperature, as a function of varying current.

functions:

Voltage_nodes : calculate the voltage at each node in the network.

find_currents : calculate the currents in each resistor and its heat production.

heat_transfer_spectral : disperses heat in the system through diffusion.

fft2d : diffusion of heat.

update_resistors_SC :  calculates the states and resistance of the resistors following the updated temperature.



