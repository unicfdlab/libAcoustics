# Contents

1. [About libAcoustics library](#About-libAcoustics-library)
2. [Available versions of the library](#Available-versions-of-the-library)
3. [Tutorials and examples](#Tutorials-and-examples)
4. [Published papers related to libAcoustics technology](#Published-materials-related-to-libAcoustics-library)
6. [For citation](#For-citation)

# About libAcoustics library
[To the contents](#Contents)

The libAcoustics - library for far-field noise computation. At current time contains 3 kind of far-field prediction models:

1. **Curle analogy**

   Based on equation 3.3 from [Curle N. (1955)](https://royalsocietypublishing.org/doi/10.1098/rspa.1955.0191)

2. **Ffowcs Williams Hawkings analogy**

     
   * **Farassat 1A integral formulation** [Brentner K. S., & Farassat F. (1998).](https://doi.org/10.2514/2.558)
   
   * **Wind-tunnel configuration: formulation GT** [Guillaume Brès et al. (2010)](https://doi.org/10.2514/6.2010-3711)

3. **CFD-BEM coupling** 

# Tutorials and examples
[To the contents](#Contents)

The libAcoustics source code tutorial, which was prepared in [Chalmers University](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2021/DebarsheeGhosh/OSCFD21_FinalPresentation.pdf)

Some examples are presented in the conference paper [Epikhin A.(2021)](https://doi.org/10.5281/zenodo.5906668)

libAcoustics library Wiki (in progress)

# Available versions of the library 
[To the contents](#Contents)

We currently support the library for the OpenFOAM+ versions. 
[Current Release](https://github.com/unicfdlab/libAcoustics/releases/tag/digitef-dev-2112)

The library is available for next OpenFOAM versions:

* OpenFOAM 3.0, 4.1 versions of the library are stored in the [openfoam-v4.1 branch](https://github.com/unicfdlab/libAcoustics/tree/openfoam-v4.1)
* OpenFOAM+ v.1812 version of the library is in the [digitef-dev-1812 branch](https://github.com/unicfdlab/libAcoustics/tree/digitef-dev-1812) of the repository
* OpenFOAM+ v.1912 version of the library is in the [digitef-dev-1912 branch](https://github.com/unicfdlab/libAcoustics/tree/digitef-dev-1912) of the repository
* OpenFOAM+ v.2012 version of the library is in the [digitef-dev-2012 branch](https://github.com/unicfdlab/libAcoustics/tree/digitef-dev-2012) of the repository
* OpenFOAM+ v.2112 version of the library is in the [digitef-dev-2112 branch](https://github.com/unicfdlab/libAcoustics/tree/digitef-dev-2112) of the repository
* OpenFOAM+ v.2212 version of the library is in the [v2212 branch](https://github.com/unicfdlab/libAcoustics/tree/v2212) of the repository
* OpenFOAM+ v.2312 version of the library is in the [v2312 branch](https://github.com/unicfdlab/libAcoustics/tree/v2312) of the repository
* OpenFOAM+ v.2412 version of the library is in the [v2412 branch](https://github.com/unicfdlab/libAcoustics/tree/v2412) of the repository

Note: Please use the v2112 or later versions because we fixed some important bugs in them that had been present in earlier implementations.


# Published materials related to libAcoustics library
[To the contents](#Contents)

## <p align="center"> >>>>> 2025 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Cavitating wake dynamics and hydroacoustics performance of marine propeller with a nozzle](https://pubs.aip.org/aip/pof/article-abstract/37/1/015128/3329276/Cavitating-wake-dynamics-and-hydroacoustics?redirectedFrom=fulltext)| --- |
|[Numerical Investigation of Unsteady Flow Characteristics of Planar Plug Nozzle in Overexpanded Conditions](https://asmedigitalcollection.asme.org/fluidsengineering/article-abstract/147/5/051203/1208662/Numerical-Investigation-of-Unsteady-Flow?redirectedFrom=fulltext)|![Jets Schlieren](https://github.com/mkraposhin/libAcoustics/blob/master/figures/2024-09-09_08-39.png)|
|[Validation of a Hybrid Approach for Wind Turbine Noise Prediction Using the Linearized Euler Equations](https://link.springer.com/chapter/10.1007/978-3-031-71926-4_16)|![The numerical problem setup](https://github.com/mkraposhin/libAcoustics/blob/master/figures/42405_2020_340_Fig1_HTML.png)|

## <p align="center"> >>>>> 2024 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Aerodynamic noise generated in three-dimensional lock-in and galloping behavior of square cylinder in high Reynolds number turbulent flows]()| --- |
|[Modeling of hydroacoustic noise from marine propellers with tip vortex cavitation]()| --- |
|[High fidelity CFD methods]()| --- |
|[Aerodynamic noise and its passive control of subcritical flow past a cylinder]()| --- |
|[Acoustic Characteristics of Laminar-Turbulent Transition in a Flat Plate Based on PSE and FW-H]()|---|
|[Low Noise Propeller Design Technology for Urban Air Mobility(UAM)](https://www.dbpia.co.kr/Journal/articleDetail?nodeId=NODE11793279)| --- |
|[Evaluation of different FW-H surfaces and modal decomposition techniques for the acoustic analysis of UAV propellers through detached eddy simulations](https://www.sciencedirect.com/science/article/pii/S1270963824000890): **Article**|![FW-H surfaces: cylinders (in red) and spheres (in blue)](https://github.com/mkraposhin/libAcoustics/blob/master/figures/1-s2.0-S1270963824000890-gr007.jpg)|
|[Hydroacoustics of turbulent flow over superhydrophobic and oscillating cylinders](https://www.sciencedirect.com/science/article/abs/pii/S002980182400088X): **Article**|![Q-criterion and flow structures](https://github.com/mkraposhin/libAcoustics/blob/master/figures/1-s2.0-S002980182400088X-gr5.jpg)|
|[Low-Noise Airfoils for Turbomachinery Applications: Two Examples of Optimization](https://www.mdpi.com/2504-186X/9/1/9): **Article**|![Isosurface for the optimal high-Reynolds airfoil](https://github.com/mkraposhin/libAcoustics/blob/master/figures/ijtpp-09-00009-g006.png)|
|[Investigation of noise reduction effects of individual and combined geometrical modifications on UAM propeller blades](https://iopscience.iop.org/article/10.1088/1742-6596/2716/1/012060/meta): **Article**|![UAM propeller blade geometry](https://github.com/mkraposhin/libAcoustics/blob/master/figures/uam_propeller_blade.png)|

## <p align="center"> >>>>> 2023 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Numerical investigation of noise suppression and amplification in forced oscillations of single and tandem cylinders in high Reynolds number turbulent flows](https://www.sciencedirect.com/science/article/abs/pii/S0307904X23000069): **Article**|![Flow fields around cylinders](https://github.com/mkraposhin/libAcoustics/blob/master/figures/1-s2.0-S0307904X23000069-gr4.jpg)|
|[AeroSPEED: a high order acoustic solver for aeroacoustic applications](https://link.springer.com/chapter/10.1007/978-3-031-40864-9_3): **Article**|![AeroSPEED](https://github.com/mkraposhin/libAcoustics/blob/master/figures/aerospeed.png)|
|[Numerical investigation of splitter plates as noise reduction techniques for tandem cylinders](https://www.researchgate.net/publication/377063820_Numerical_investigation_of_splitter_plates_as_noise_reduction_techniques_for_tandem_cylinders): **Article**|![The mesh around cylinders](https://github.com/mkraposhin/libAcoustics/blob/master/figures/two_cyls_mesh.png)|
|[Comparison of Porous and Direct Volumetric Integration FW-H Formulation for Acoustic Prediction](https://onepetro.org/ISOPEIOPEC/proceedings-abstract/ISOPE23/All-ISOPE23/524676): **Article**| --- |
|[Prediction & Active Control of Multi-Rotor Noise](https://commons.erau.edu/edt/732/): **MSc Thesis**|![Schematic of various multi-rotor noise sources](https://github.com/mkraposhin/libAcoustics/blob/master/figures/noise-sources.png)|
|[Noise prediction of a single-stream under-expanded jet in OpenFOAM](https://arc.aiaa.org/doi/abs/10.2514/6.2023-3349): **Article**|---|
|[Numerical Simulation of Supersonic Jet Noise Using Open Source Software](https://link.springer.com/chapter/10.1007/978-3-031-36030-5_24): **Article**| --- |
|[Evaluation of Synthetic Jet Flow Control Technique for Modulating Turbulent Jet Noise](https://www.mdpi.com/2311-5521/8/4/110): **Article**|![Coherent structures obtained using Q-criterion](https://github.com/mkraposhin/libAcoustics/blob/master/figures/fluids-08-00110-g013.png)|
|[Numerical Investigation on Aerodynamic Noise of Flow past a Cylinder with Different Spanwise Lengths](https://aip.scitation.org/doi/abs/10.1063/5.0139731):  **Article**|![A Q-criterion for different size cylinders](https://github.com/mkraposhin/libAcoustics/blob/master/figures/Guanjiang_Chen_1.png)|
|[Aerodynamic and Acoustic Analysis of Multiple Rotors in Different Configurations](https://www.researchgate.net/publication/367313189_Aerodynamic_and_Acoustic_Analysis_of_Multiple_Rotors_in_Different_Configurations): **Article**|---|
|[CFD-based Aerodynamic and Aeroacoustic Analysis of Large Payload Multi-Copter Rotors](https://www.researchgate.net/publication/367317951_CFD-based_Aerodynamic_and_Aeroacoustic_Analysis_of_Large_Payload_Multi-Copter_Rotors): **Article**|![TKE plots and spectra](https://github.com/mkraposhin/libAcoustics/blob/master/figures/Naina_Pisharoti.png)|
|[Simulations of multi-rotor interaction noise at hovering & forward flight conditions](https://www.researchgate.net/publication/367322923_Simulations_of_multi-rotor_interaction_noise_at_hovering_forward_flight_conditions): **Article**|---|
|[Combined propeller blade geometrical modifications to reduce noise emissions of urban air mobility vehicles](https://repositum.tuwien.at/handle/20.500.12708/188143)| --- |

## <p align="center"> >>>>> 2022 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Aeroacoustic characterization of a 3D organ pipe](https://www.politesi.polimi.it/handle/10589/197383)| --- |
|[Simetrik Kanat Profilinin Aerodinamik ve Aeroakustik Olarak İncelenmesi](https://www.researchgate.net/profile/Mehmet-Alim-2/publication/379825150_Simetrik_Kanat_Profilinin_Aerodinamik_ve_Aeroakustik_Olarak_Incelenmesi/links/661d174866ba7e2359dcf236/Simetrik-Kanat-Profilinin-Aerodinamik-ve-Aeroakustik-Olarak-Incelenmesi.pdf)| --- |
|[Towards a low-noise axial fan for automotive applications](https://www.researchgate.net/publication/366519846_Towards_a_low-noise_axial_fan_for_automotive_applications): **Article**|![Test rig isometry and sketches](https://github.com/mkraposhin/libAcoustics/blob/master/figures/blade-test-rig.png)|
|[Numerical Study of Cavitation Noise Around NACA66 (MOD) Hydrofoil with Direct Volume Integration](https://onepetro.org/IJOPE/article-abstract/33/03/294/533531/Numerical-Study-of-Cavitation-Noise-Around-NACA66?redirectedFrom=fulltext): **Article**|![The cavity near hydrofoil](https://github.com/mkraposhin/libAcoustics/blob/master/figures/hydrofoil_w_cavity.png)|
|[Numerical Study of the Impact of Fluid–Structure Interaction on Flow Noise over a Rectangular Cavity](https://www.mdpi.com/1996-1073/15/21/8017): **Article**|![Modes of a cavity walls](https://github.com/mkraposhin/libAcoustics/blob/master/figures/cav_walls_modes.png)|
|[Validating Confined Flame Noise Simulation Using External Sensor](https://www.mdpi.com/1424-8220/22/20/8039): **Article**|![Contour of flame temperature](https://www.mdpi.com/sensors/sensors-22-08039/article_deploy/html/images/sensors-22-08039-g006-550.jpg)|
|[Numerical Study of Cavitation Noise Around 2-D Hydrofoil with Direct Volume Integration](https://onepetro.org/ISOPEIOPEC/proceedings-abstract/ISOPE22/All-ISOPE22/ISOPE-I-22-497/494136): **Article**| --- |
|[Prediction of the flow and acoustic fields generated by an isolated propeller](https://www.researchgate.net/publication/363298683_Prediction_of_the_flow_and_acoustic_fields_generated_by_an_isolated_propeller): **Article**|![A propeller isometric view](https://github.com/mkraposhin/libAcoustics/blob/master/figures/daSilvaetal_FIA_2020_22.png)|
|[SIMULATION OF THE FLOW FIELD AND FAR-FIELD NOISE OF AN UAV PROPELLER](https://www.researchgate.net/profile/Filipe-Dutra-Da-Silva/publication/366047502_Simulation_of_the_flow_field_and_far-field_noise_from_an_UAV_propeller/links/663f8f2335243041539ad217/Simulation-of-the-flow-field-and-far-field-noise-from-an-UAV-propeller.pdf)| --- |
|[VALIDATION AND VERIFICATION OF THE ACOUSTIC MODEL ON THE PROBLEM OF FLOW AROUND A ROUND CYLINDER](https://www.lki.smtu.ru/en/article/68/)| --- |
|[Simulation of Combustion Noise of Premixed Flames in OpenFOAM](https://search.proquest.com/openview/bc825d9e785d2a6e5c052d97e9f07f19/1?pq-origsite=gscholar&cbl=18750&diss=y)| --- |

## <p align="center"> >>>>> 2021 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Wall-modeled Large-Eddy Simulations for Trailing-Edge Turbulent Boundary Layer Noise Prediction
](https://rave.ohiolink.edu/etdc/view?acc_num=osu1610043130955909)| --- |
|[Evaluación de la técnica de control de flujo tipo chorro sintético para la modulación de niveles de ruido en flujos de inyección turbulenta](https://repositorio.unal.edu.co/handle/unal/81180)| --- |
|[Quadcopter drones swarm aeroacoustics](https://doi.org/10.1063/5.0052505)|![Flow fields around quadcopters](https://github.com/mkraposhin/libAcoustics/blob/master/figures/images-qcs.png)|
|[Computational aeroacoustics of quadcopter drones](https://doi.org/10.1016/j.apacoust.2022.108738)|![Flow fields](https://github.com/mkraposhin/libAcoustics/blob/master/figures/1-s2.0-S0003682X22001128-gr5.jpg)|
|[Implementation of the FWH aero-acoustic analogy for sector analysis of an axi-symmetric turbomachine](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2021/DebarsheeGhosh/OSCFD21_FinalPresentation.pdf): **Tutorial**|![wedge domain](https://github.com/unicfdlab/libAcoustics/blob/master/figures/fan-wedge-domain.png)|
|[Prediction of infrasound emission from horizontal axus wind turbines](https://github.com/unicfdlab/PhDTheses/blob/main/Thesis_Copy_edited-1.pdf): **MSc Thesis**|![Computational method scketch](https://github.com/unicfdlab/libAcoustics/blob/master/figures/Infra_MScThesis.png)|
|[Influence of the elastic cavity walls on cavity flow noise](https://www.researchgate.net/publication/356458463_Influence_of_the_elastic_cavity_walls_on_cavity_flow_noise):  **Article** |![flexible cavity instant velocity field](https://github.com/unicfdlab/libAcoustics/blob/master/figures/flex-cavity-U.png)|
|[A Computational Study on the Aeroacoustics of a Multi-Rotor Unmanned Aerial System](https://www.researchgate.net/publication/355396901_A_Computational_Study_on_the_Aeroacoustics_of_a_Multi-Rotor_Unmanned_Aerial_System):  **Article** |![Unmanned Aerial System](https://github.com/unicfdlab/libAcoustics/blob/master/figures/Unmanned-Aerial-System.png)|
|[Wall-modeled large-eddy simulation of a trailing-edge serration–finlet configuration](https://www.researchgate.net/publication/352381522_Wall-Modeled_Large-Eddy_Simulation_of_a_Trailing-Edge_Serration-Finlet_Configuration):  **Article** |![Streamwise vector plots in the serration–finlet configuration](https://aip.scitation.org/na101/home/literatum/publisher/aip/journals/content/adv/2021/adv.2021.11.issue-6/5.0049181/20210613/images/medium/5.0049181.figures.online.f30.jpg)|
|[The Eulerian–Lagrangian Approach for the Numerical Investigation of an Acoustic Field Generated by a High-Speed Gas-Droplet Flow](https://www.mdpi.com/2311-5521/6/8/274):  **Article** | ![Jet with particles Logo](https://www.mdpi.com/fluids/fluids-06-00274/article_deploy/html/images/fluids-06-00274-ag-550.jpg)|
|[Validation of the developed open-source library for far-field noise prediction](https://www.researchgate.net/publication/354447445_Validation_of_the_developed_open-source_library_for_far-field_noise_prediction):  **Article** |![Domain parameters: (a) boundary conditions; (b) mesh configuration and FW-H surface.](https://github.com/unicfdlab/libAcoustics/blob/master/figures/Validation-libAcoustics.jpg)|
|[A proposed wavy shield for suppression of supersonic jet noise utilizing reflections](https://www.researchgate.net/publication/348451623_A_proposed_wavy_shield_for_suppression_of_supersonic_jet_noise_utilizing_reflections)|![Acoustic fields](https://github.com/unicfdlab/libAcoustics/blob/master/figures/Capture.PNG)|
|[Numerical Study of Flow-Induced Noise around Cylinder and Rod-NACA0012 Hydrofoil](https://onepetro.org/ISOPEIOPEC/proceedings-abstract/ISOPE21/All-ISOPE21/ISOPE-I-21-1191/464451)|---|
|[Simulations of Noise Generated by Rotor-Rotor Interactions at Static Conditions](https://arc.aiaa.org/doi/abs/10.2514/6.2021-1986)|---|
|[Investigating the Impact of Water Injection on Noise Generation During Rocket Lift-Off](https://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1598708&dswid=4558)|---|
|[Análise e validação de modelo de simulação numérica para previsão de ruído aerodinâmico utilizando OpenFOAM e Libacoustics](https://repositorio.ufsc.br/handle/123456789/228458?show=full): **MSc Thesis**|![acoustic field around tandem](https://github.com/unicfdlab/libAcoustics/blob/master/figures/Tandem-acoustic-field.png)|
|[PREDICTION OF AERODYNAMIC NOISE GENERATED BY CYLINDERS IN UNIFORM FLOW USING OPENFOAM](https://www.researchgate.net/profile/Filipe-Dutra-Da-Silva/publication/356642546_Prediction_of_aerodynamic_noise_generated_by_cylinders_in_uniform_flow_using_OpenFOAM/links/663f904f7091b94e931ddae8/Prediction-of-aerodynamic-noise-generated-by-cylinders-in-uniform-flow-using-OpenFOAM.pdf)| --- |
|[Numerical simulation of rotating systems in fluid dynamics using free access computational tools: Application to Balance Ring](https://somim.org.mx/memorias/memorias2021/articulos/A5_55.pdf)| --- |
|[Prediction of Infrasound Emission from Horizontal Axis Wind Turbines](https://search.proquest.com/openview/0d10d5ddce2a16a3fa6794babea9dcec/1?pq-origsite=gscholar&cbl=18750&diss=y)| --- |
|[Modeling of a ventilated cavity behind a streamlined body](https://www.researchgate.net/profile/Nf-Dimitrieva/publication/354611649_MODELING_OF_A_VENTILATED_CAVITY_BEHIND_A_STREAMLINED_BODY/links/61ffad47870587329e9720f7/MODELING-OF-A-VENTILATED-CAVITY-BEHIND-A-STREAMLINED-BODY.pdf)| --- |

## <p align="center"> >>>>> 2020 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Numerical Study of Iced Airfoil Aeroacoustics Using IDDES](https://www.researchgate.net/publication/342020749_Numerical_Study_of_Iced_Airfoil_Aeroacoustics_Using_IDDES):  **Article** | ![Q-criterion around wing](https://github.com/unicfdlab/libAcoustics/blob/master/figures/L0.2Pc_Qcriterion.png) |
|[Jet Noise in Airframe Integration and Shielding](https://www.mdpi.com/2076-3417/10/2/511):  **Article** |![Emittance of acoustic waves by jet near wing](https://www.mdpi.com/applsci/applsci-10-00511/article_deploy/html/images/applsci-10-00511-g003-550.jpg)|
|[Coupled CFD-based Shape Optimization of a Wing of Reusable Space Vehicle of Tourist Class](https://www.researchgate.net/publication/338372884_Coupled_CFD-based_Shape_Optimization_of_a_Wing_of_Reusable_Space_Vehicle_of_Tourist_Class/figures):  **Article**|![Computational scheme](https://github.com/unicfdlab/libAcoustics/blob/master/figures/wingW640.jpg)|
|[Transient cavitating flow structure and acoustic analysis of a hydrofoil with whalelike wavy leading edge](https://www.researchgate.net/publication/341124148_Transient_cavitating_flow_structure_and_acoustic_analysis_of_a_hydrofoil_with_whalelike_wavy_leading_edge):  **Article**|![Flow fields around wing](https://github.com/unicfdlab/libAcoustics/blob/master/figures/hydrofoilW640.jpg)|
|[Prediction of the Free Jet Noise Using Quasi-gas Dynamic Equations and Acoustic Analogy](https://link.springer.com/chapter/10.1007/978-3-030-50436-6_16): **Article**|![QGDFoam instant jet velocities](https://media.springernature.com/lw785/springer-static/image/chp%3A10.1007%2F978-3-030-50436-6_16/MediaObjects/500810_1_En_16_Fig5_HTML.png)|
|[Simulations of Broadband Noise of a Small UAV Propeller](https://arc.aiaa.org/doi/abs/10.2514/6.2020-1493): **Article**|---|
|[High-Fidelity Simulations of Noise Generation in a Propeller-Driven Unmanned Aerial Vehicle](https://arc.aiaa.org/doi/abs/10.2514/1.J059117): **Article**|---|
|[Finite Element Modelling of a Flow-Acoustic Coupling in Unbounded Domains](https://journals.pan.pl/Content/118240/PDF/aoa.2020.135251.pdf)| --- |


## <p align="center"> >>>>> 2019 <<<<<  </p>

| Title | Description |
|------|-------------|
|[The numerical simulation of compressible jet at low Reynolds number using OpenFOAM](https://www.researchgate.net/publication/337116883_The_numerical_simulation_of_compressible_jet_at_low_Reynolds_number_using_OpenFOAM):  **Article**|![Sound pressure level directivity distributions for M = 0.9, Re = 3600](https://github.com/unicfdlab/libAcoustics/blob/master/figures/easLibAcoustics.jpg)|
|[Wing airfoil selection and optimization for the tourist class reusable space vehicle](https://www.researchgate.net/publication/335093847_Wing_airfoil_selection_and_optimization_for_the_tourist_class_reusable_space_vehicle/citations):  **Article**|---|
|[Analysis Methods and Design Measures for the Reduction of Noise and Vibration Induced by Marine Propellers](http://pub.dega-akustik.de/ICA2019/data/articles/001556.pdf): **Article**|![Vortices after propeller thruster](https://github.com/unicfdlab/libAcoustics/blob/master/figures/PT-vortices.png)|
|[Prediction of Noise Associated with an Isolated UAV Propeller](https://commons.erau.edu/edt/463/): **MSc Thesis**|![Velocity vectors around wing](https://github.com/unicfdlab/libAcoustics/blob/master/figures/wing-velocityfield.png)|

## <p align="center"> >>>>> 2018 <<<<<  </p>
| Title | Description |
|------|-------------|
|[Multiscale approach for simulation of complex transient processes of fluid flows in technical systems](https://ispranproceedings.elpub.ru/jour/issue/download/84/55#page=291)| --- |

## <p align="center"> >>>>> 2016 <<<<<  </p>

| Title | Description |
|------|-------------|
|[Broadband noise prediction using large eddy simulation and a frequency domain method](https://www.researchgate.net/publication/309877358_Broadband_noise_prediction_using_large_eddy_simulation_and_a_frequency_domain_method):  **Article**|![Span-wise vorticity contour](https://ars.els-cdn.com/content/image/1-s2.0-S0003682X16304145-gr11.jpg)|


# For citation
[To the contents](#Contents)

If you have found the software useful for your research, please cite next sources:

* Epikhin, A., Evdokimov, I., Kraposhin, M., Kalugin, M., Strijhak, S. Development of a Dynamic Library for Computational Aeroacoustics Applications Using the OpenFOAM Open Source Package // Procedia Computer ScienceVolume 66, 2015, Pages 150-157
https://www.sciencedirect.com/science/article/pii/S1877050915033670 , DOI: 10.1016/j.procs.2015.11.018
* A. Epikhin, Validation of the developed open-source library for far-field noise prediction, in: Proceedings of the 7th International Congress on Sound and Vibration, 2021. DOI: 10.5281/zenodo.5906668
* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3878439.svg)](https://doi.org/10.5281/zenodo.3878439) 
