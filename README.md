# Active clock model and active XY model

Codes used in the scientific publication:<br>
S. Chatterjee, M. Mangeat, and H. Rieger, <a href='https://iopscience.iop.org/article/10.1209/0295-5075/ac6e4b'>Polar flocks with discretized directions: the active clock model approaching the Vicsek model</a>, EPL <b>138</b>, 41001 (2022). Preprint available on <a href='https://arxiv.org/abs/2203.01181'>arXiv</a>.

## q-state active clock model

C++ code on numerical simulation of the q-state active clock model.<br>
Exportations: dynamics of flocking; time evolution of the total magnetization and the mean-square displacement; number and magnetization fluctuations.<br>
Compile: g++ active_clock_model.cpp -lgsl -lgslcblas -lm -O3 -s -o active_clock_model.out.<br>
Run: ./active_clock_model.out -parameter=value.<br>
List of parameters: q, beta, epsilon, rho0, LX, LY, init, RAN, tmax (details as comments in the code).

## active XY model

C++ code on numerical simulation of the active XY model.<br>
Exportations: dynamics of flocking; time evolution of the total magnetization and the mean-square displacement; number and magnetization fluctuations.<br>
Compile: g++ active_XY_model.cpp -lgsl -lgslcblas -lm -O3 -s -o active_XY_model.out.<br>
Run: ./active_XY_model.out -parameter=value.<br>
List of parameters: beta, epsilon, rho0, LX, LY, init, RAN, tmax (details as comments in the code).
