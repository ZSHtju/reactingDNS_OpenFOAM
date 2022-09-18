# 1DNHPTNH3Simulations
> This is a tutorial for one dimentional NHPT-NH3 flame simulations. It has three base cases in the `CASES` folders and a set of direct numerical simulation (DNS) code in the `CASES` folder. This instruction shows how to use the code and run the designed simulations planned in the `NHPTNH3BoundaryConditions.xlsx` file.

## Prerequisite
- OpenFOAM-7

## Compiling codes
This compiling code process is required before runing the cases, and you only need to do it once.

Download the code and unzip it. Typing the following lines to compile it.
```
cd 1DNHPTNH3Simulations/CODES/CCMChemistryModel
wmake
cd ../reactingDNS/solvers/reactingDNS
wmake
```

## Runing cases
The `NHPTNH3BoundaryConditions.xlsx` file shows the desired 1D simulations. It has 6 sheets, baseCase, EQ0.0, EQ0.2, ..., EQ1.0, showing 6 different premixed NH3/air equivalence ratios (EQ). Each sheet has 70 lines shows the engine in-cylinder conditions ($T$, $P$, and species mass fractions $Y$) at different engine crank angles. Each line correspond to an 1D case, where the boundary conditions ($T$, $P$, and mass fractions of NH3, O2, N2 in the 0 folder) are set for the premixed NH3/air mixture.

To run a single dual-fuel 1D case, you need to follow two step, i.e., initialization and formal simulations. The initialization step is used to generate an initial distribution of the non-reacting C7H16/NH3 mixing part. The formal simulations are used to simulate the ignition, flame propogation and autoignition of a typical dual fuel combustion under constant volume and constant pressure conditions.

### Initialization
In this step, a laminar counter flow case is configured and performed to generate temperature (T) and species mass fractions (Y) profiles, cf. the folowing profiles.

![Profiles](CASES/1.Initialization/counterFlowSimpleMechanism/mixingLayerProfiles.png)

The case in located in `CASES/1.Initialization/counterFlowSimpleMechanism` and it has an 1D domain with 1000 computational cells and a dimension of H = 0.01 m in $x$-direction. The left side is fuel inlet and the right side is air inlet.  It is a non-reacting mixing with a strain rate $a$ = $U$/H, in which $U$ is the inlet velocity and is set to be 0.1 m/s (see `0/U` file in the `CASES/1.Initialization/counterFlowSimpleMechanism` folder). Thus the strain rate we studied is $a$ = 10.

When you setup a new Initialization case, you need to copy the `counterFlowSimpleMechanism` case in the `CASES/1.Initialization/counterFlowSimpleMechanism` folder and change the following lines according to the `NHPTNH3BoundaryConditions.xlsx` file.
1. In `0/p` file, line 20;
2. In `0/T` file, line 20, 32, 37 and 38;
3. In `0/NH3`, `0/C7H16`, `0/O2`, and `0/N2` files, line 32.

### FormalSimulations
There are two similar cases in the `CASES/2.FormalSimulations` folder, i.e., `1D_NH3_NHPT_ConstP` and `1D_NH3_NHPT_ConstV`, namely a constant pressure and a constant volume configuration, respectively. The only difference is in the `0/T` and `0/p` files. In the `1D_NH3_NHPT_ConstP` case, `0/T` and `0/p` files are set to be an open domain while they are set to be a closed domain in the the `1D_NH3_NHPT_ConstV` . Thus, the `1D_NH3_NHPT_ConstP` case will have a constant pressure during the simulation and the `1D_NH3_NHPT_ConstV` case will have a increasing pressure during the simulation.

When you setup a new FormalSimulations case, you need to copy the `1D_NH3_NHPT_ConstP` and `1D_NH3_NHPT_ConstV` cases. Let us start from `1D_NH3_NHPT_ConstP` first, you need to change the following lines according to the `NHPTNH3BoundaryConditions.xlsx` file.
1. In `0/p` file, line 20, 36, 42, 48 and 54;
2. In `0/T`, `0/NH3`, `0/C7H16`, `0/O2`, and `0/N2` files, line 23-1022 (the first 1000 cells, copy from the Initialization case), 1023-3022 (the last 2000 cells, which are premixed mixtures, use the same value as in `NHPTNH3BoundaryConditions.xlsx`);

For the `1D_NH3_NHPT_ConstV` case, you need to change the following lines
1. In `0/p` file, line 20;
2. In `0/T`, `0/NH3`, `0/C7H16`, `0/O2`, and `0/N2` files, line 23-1022 (the first 1000 cells, copy from the Initialization case), 1023-3022 (the last 2000 cells, which are premixed mixtures, use the same value as in `NHPTNH3BoundaryConditions.xlsx`);
