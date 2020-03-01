# Overhead_Power_Lines
Simulating overhead power lines in phase domain
Bachelor thesis in North China Electric Power University
Completed in June 19 2018

Abstract

Overhead line is an important power transmission equipment in power system. It is necessary to establish the broadband model of overhead lines in order to improve the precision in harmonic and transient simulation. On the whole, the frequency-dependent model of power transmission line used for dynamic characteristic simulation of power grid is based on the frequency-dependent parameters. A rational function is used to fit the frequent characteristics of the parameters by the method of Vector Fitting. Then the equivalent circuit is formed and finally used in the time domain. J Marti’s FD-Line model and its improved model are representative model applied to the electromagnetic transient simulation software. The existing ADPSS (Advanced Digital Power System Simulator) is still deficient in some aspects, such as the modeling of asymmetric multiphase transmission lines, the frequent characters of phase-mode transformation matrix and the direct modeling in phase domain. Therefore, it is important to study the phase domain broadband modeling of overhead line to satisfy the demand of electromagnetic transient simulation based on the theory of electromagnetic field and circuit analysis.

The paper is based on the ULM (Universal Line Model), which is different from J Marti’s FD-Line model. The ULM is modeled directly in the phase domain, avoiding the error caused by the constant real phase-mode transformation matrix and making the model more accurate in the frequency range of interest and the results more precise in time domain. The paper’s work is as follows: lots of literature are investigated firstly. Then the details of the model are studied and a detailed modeling process and algorithm flow are given. Lastly, an example is given to verify the effectiveness of the model.
