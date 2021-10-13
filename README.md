# Contents

1. [libAcoustics library brief](#libAcoustics-library-brief)
2. [Published papers related to libAcoustics technology](#Published-papers-related-to-libAcoustics-technology)
3. [Available versions of the library](#Available-versions-of-the-library)
4. [For citation](#For-citation)

# libAcoustics library brief
[To the contents](#Contents)

libAcoustics - library for far-field noise computation. At current time contains 3 kind of far-field predictions:

1. **Curle analogy**

2. **Ffowcs Williams Hawkings analogy**

3. **CFD-BEM coupling** 

# Available versions of the library 
[To the contents](#Contents)

The library is available for next OpenFOAM versions:

* OpenFOAM 3.0, 4.1 versions are stored in the master branch
* OpenFOAM+ v.1812 is in the digitef-dev-1812 branch of the repository
* OpenFOAM+ v.1912 is in the digitef-dev-1912 branch of the repository
* OpenFOAM+ v.2012 is in the digitef-dev-2012 branch of the repository

# Published papers related to libAcoustics library
[To the contents](#Contents)

| Title | Description |
|------|-------------|
| [Numerical Study of Iced Airfoil Aeroacoustics Using IDDES](https://www.researchgate.net/publication/342020749_Numerical_Study_of_Iced_Airfoil_Aeroacoustics_Using_IDDES):  **Article** | ![Q-criterion around wing](https://github.com/unicfdlab/libAcoustics/blob/master/L0.2Pc_Qcriterion.png) |
|[Wall-modeled large-eddy simulation of a trailing-edge serration–finlet configuration](https://www.researchgate.net/publication/352381522_Wall-Modeled_Large-Eddy_Simulation_of_a_Trailing-Edge_Serration-Finlet_Configuration):  **Article** |![Streamwise vector plots in the serration–finlet configuration](https://aip.scitation.org/na101/home/literatum/publisher/aip/journals/content/adv/2021/adv.2021.11.issue-6/5.0049181/20210613/images/medium/5.0049181.figures.online.f30.jpg)|
|[The Eulerian–Lagrangian Approach for the Numerical Investigation of an Acoustic Field Generated by a High-Speed Gas-Droplet Flow](https://www.mdpi.com/2311-5521/6/8/274):  **Article** | ![Jet with particles Logo](https://www.mdpi.com/fluids/fluids-06-00274/article_deploy/html/images/fluids-06-00274-ag-550.jpg)|
|[Jet Noise in Airframe Integration and Shielding](https://www.mdpi.com/2076-3417/10/2/511):  **Article** |![Emittance of acoustic waves by jet near wing](https://www.mdpi.com/applsci/applsci-10-00511/article_deploy/html/images/applsci-10-00511-g003-550.jpg)|
|[Validation of the developed open-source library for far-field noise prediction](https://www.researchgate.net/publication/354447445_Validation_of_the_developed_open-source_library_for_far-field_noise_prediction):  **Article** |![Domain parameters: (a) boundary conditions; (b) mesh configuration and FW-H surface.](https://www.researchgate.net/profile/Andrey-Epikhin/publication/354447445/figure/fig1/AS:1065774030548994@1631111738797/Domain-parameters-a-boundary-conditions-b-mesh-configuration-and-FW-H-surface_W640.jpg)|
|[Broadband noise prediction using large eddy simulation and a frequency domain method](https://www.researchgate.net/publication/309877358_Broadband_noise_prediction_using_large_eddy_simulation_and_a_frequency_domain_method)):  **Article**|![Span-wise vorticity contour](https://ars.els-cdn.com/content/image/1-s2.0-S0003682X16304145-gr11.jpg)|
|[The numerical simulation of compressible jet at low Reynolds number using OpenFOAM](https://www.researchgate.net/publication/337116883_The_numerical_simulation_of_compressible_jet_at_low_Reynolds_number_using_OpenFOAM):  **Article**|![Sound pressure level directivity distributions for M = 0.9, Re = 3600](https://www.researchgate.net/publication/337116883/figure/fig2/AS:823141379092492@1573263606747/Sound-pressure-level-directivity-distributions-for-M-09-Re-3600-and-PCF-solver-with_W640.jpg)|
|[Wing airfoil selection and optimization for the tourist class reusable space vehicle](https://www.researchgate.net/publication/335093847_Wing_airfoil_selection_and_optimization_for_the_tourist_class_reusable_space_vehicle/citations):  **Article**|---|
|[Coupled CFD-based Shape Optimization of a Wing of Reusable Space Vehicle of Tourist Class](https://www.researchgate.net/publication/338372884_Coupled_CFD-based_Shape_Optimization_of_a_Wing_of_Reusable_Space_Vehicle_of_Tourist_Class/figures):  **Article**|![Computational scheme](https://www.researchgate.net/publication/338372884/figure/fig1/AS:843306078568450@1578071245982/Design-variables-of-the-wing_W640.jpg)|
|[Transient cavitating flow structure and acoustic analysis of a hydrofoil with whalelike wavy leading edge](https://www.researchgate.net/publication/341124148_Transient_cavitating_flow_structure_and_acoustic_analysis_of_a_hydrofoil_with_whalelike_wavy_leading_edge):  **Article**|![Flow fields around wing](https://www.researchgate.net/profile/Zehao-Li-3/publication/341124148/figure/fig12/AS:890182240837633@1589247393515/Pressure-coefficient-C-p-and-limiting-streamline-colored-by-velocity-in-flow_W640.jpg)|
|[A proposed wavy shield for suppression of supersonic jet noise utilizing reflections](https://www.researchgate.net/publication/348451623_A_proposed_wavy_shield_for_suppression_of_supersonic_jet_noise_utilizing_reflections)|![Acoustic fields](https://github.com/unicfdlab/libAcoustics/blob/master/Capture.PNG)|
|[Numerical Study of Flow-Induced Noise around Cylinder and Rod-NACA0012 Hydrofoil](https://onepetro.org/ISOPEIOPEC/proceedings-abstract/ISOPE21/All-ISOPE21/ISOPE-I-21-1191/464451)|---|
|[Simulations of Noise Generated by Rotor-Rotor Interactions at Static Conditions](https://arc.aiaa.org/doi/abs/10.2514/6.2021-1986)|---|
|[Investigating the Impact of Water Injection on Noise Generation During Rocket Lift-Off](https://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1598708&dswid=4558)|---|
|[Prediction of the Free Jet Noise Using Quasi-gas Dynamic Equations and Acoustic Analogy](https://link.springer.com/chapter/10.1007/978-3-030-50436-6_16): **Article**|![QGDFoam instant jet velocities](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-50436-6_16/MediaObjects/500810_1_En_16_Fig5_HTML.png)|
|[Simulations of Broadband Noise of a Small UAV Propeller](https://arc.aiaa.org/doi/abs/10.2514/6.2020-1493): **Article**|---|
|[Analysis Methods and Design Measures for the Reduction of Noise and Vibration Induced by Marine Propellers](http://pub.dega-akustik.de/ICA2019/data/articles/001556.pdf): **Article**|![Vortices after propeller thruster](https://github.com/unicfdlab/libAcoustics/blob/master/PT-vortices.png)|
|[Prediction of Noise Associated with an Isolated UAV Propeller](https://commons.erau.edu/edt/463/): **MSc Thesis**|![Velocity vectors around wing](https://github.com/unicfdlab/libAcoustics/blob/master/wing-velocityfield.png)|
|[High-Fidelity Simulations of Noise Generation in a Propeller-Driven Unmanned Aerial Vehicle](https://arc.aiaa.org/doi/abs/10.2514/1.J059117): **Article**|---|
|[Análise e validação de modelo de simulação numérica para previsão de ruído aerodinâmico utilizando OpenFOAM e Libacoustics](https://repositorio.ufsc.br/handle/123456789/228458?show=full): **MSc Thesis**|![acoustic field around tandem](https://github.com/unicfdlab/libAcoustics/blob/master/Tandem-acoustic-field.png)|



# For citation
[To the contents](#Contents)

If you have found the software useful for your research, please cite next sources:

* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3878439.svg)](https://doi.org/10.5281/zenodo.3878439) 
 
* Epikhin, A., Evdokimov, I., Kraposhin, M., Kalugin, M., Strijhak, S. Development of a Dynamic Library for Computational Aeroacoustics Applications Using the OpenFOAM Open Source Package // Procedia Computer ScienceVolume 66, 2015, Pages 150-157
https://www.sciencedirect.com/science/article/pii/S1877050915033670 , DOI: 10.1016/j.procs.2015.11.018
