# RotaryDryerDesigner: A MATLAB App for Preliminary Sizing and Axial Simulation of Direct-Heated Rotary Dryers

## Authors

Carlos Zalazar Oliva <sup>a</sup>  
Liomnis Osorio <sup>b</sup>  
Jorge Silva <sup>c</sup>  
Yoalbys Retirado-Mediaceja <sup>d</sup>  
Ever Góngora Leyva <sup>e</sup>  
Dina Álvarez Reyes <sup>a</sup>  
Deynier Montero <sup>d</sup> 
Joelmis Ramirez Oliva d</sup> 

<sup>*</sup> Corresponding authors  

## Affiliations

<sup>a</sup> Doctoral Program in Engineering, Macro-Faculty of Engineering (UFRO–UBB–UTALCA Consortium), Concepción 4051381, Chile  
<sup>b</sup> Department of Industrial Processes, Faculty of Engineering, Universidad Católica de Temuco, Rudecindo Ortega 2950, Temuco 4780000, Chile  
<sup>c</sup> Faculty of Engineering, Universidad del Bío-Bío, Concepción 4051381, Chile  
<sup>d</sup> Universidad de Moa “Dr. Antonio Núñez Jiménez”, Avenida Calixto García Iñiguez No. 15, Moa 83310, Cuba  
<sup>e</sup> Engineering and Physics Department, Faculty of Science and Engineering, Bindura University of Science Education, Bindura, Zimbabwe  


## Overview
RotaryDryerDesigner is a MATLAB-based application developed for the preliminary sizing and one-dimensional axial simulation of direct-heated rotary dryers. Through an interactive App Designer interface, the software integrates process definition, solid-property management, drying-kinetics analysis, and numerical solution of coupled heat and mass balances within a single engineering workflow. Users can enter operating conditions, fit drying models from experimental data, estimate temperature-dependent kinetic parameters, and compute the axial evolution of moisture and temperature variables along the dryer. The application generates design-oriented outputs including estimated dryer length, outlet temperatures, residence time, axial profiles, and selected performance indicators under prescribed operating conditions and target final moisture content.


## Installation
### Option 1: Run from source
1. Clone or download this repository.
2. Add the repository root to the MATLAB path.
3. Open the App Designer file:
   - `app/RotaryDryerDesigner.mlapp`

### Option 2: Executable version
If an installer package is provided, run the installer and follow the on-screen steps.

## Getting started
1. Open the application.
2. Define the process and operating conditions in the Input tab.
3. Load or define the required material properties in the Properties tab.
4. Import experimental drying data and fit the kinetic model in the Kinetics tab.
5. Enter design and numerical parameters in the Results tab.
6. Run the calculation and inspect the predicted axial profiles and summary outputs.

## Input data
The software may require:
- Process and operating conditions.
- Solid and gas properties.
- Drying-kinetics data from experiments.
- User-defined design specifications such as target outlet moisture content.

## Outputs
Depending on the case definition, the software provides:
- Estimated dryer length.
- Outlet solid temperature.
- Outlet gas temperature.
- Final solid moisture content.
- Residence time.
- Axial temperature profiles.
- Axial moisture profiles.
- Summary tables and exportable reports.

## Scope
RotaryDryerDesigner is intended for preliminary analysis of direct-heated rotary dryers using a one-dimensional axial representation. Applications involving more complex flow structures, alternative dryer configurations, or multidimensional effects may require model extensions.


## License
This repository is currently shared for academic and research purposes. Please contact the authors for permissions regarding use or distribution.

## Support
Email: rotarydryerdesigner@gmail.com