# Algorithm for Heartbeat Detection

## Getting Started
This project explains and demonstrates an afresh technique that can determine inter-beat-interval (IBI) of heartbeats and time interval of each breathing cycle caused by a human body. An aim of the project is same as ECG with some advantages like being wireless and more convenient. This extracted data can be used to predict vital signs of humans. For extraction of the heartbeats from acquired raw data an algorithm is developed and implemented in Python_2.7.

Further details can be found in report.

## Few Highlights 
* Successfully measured Heart Rate Variability and breathing rate of a person 50 cm away using Radar
* Developed an algorithm based on Joint Optimization to extract HRV from RF signals with the aim of
recognizing emotions of a person
* Developed an algorithm to identify the frequency of a moving object using Continuous Wave Radar
*	Employed Cubic Spline interpolation

## about
Textfiles contains radar's raw data of a moving plate with maximum displacement and frequency. Original heartbeat data is not included here. Each version has different updates which commented in the code itselves.

## Packages
sudo apt-get install python-numpy python-scipy python-matplotlib

## References
* Mingmin Zhao, Fadel Adib, Dina Katabi. “Emotion Recognition using Wireless Signals”
* F. Adib, H. Mao, Z. Kabelac, D. Katabi, and R. C. Miller. Smart homes that monitor breathing and heart rate. ACM, 2015
* F. Adib and D. Katabi. See through walls with wifi! In Proceedings of the ACM SIGCOMM 2013, pages 75–86, New York, NY, USA, 2013. ACM
* http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
* S. McKinley and M. Levine. Cubic spline interpolation. College of the Redwoods, 45(1):1049–1060, 1998
