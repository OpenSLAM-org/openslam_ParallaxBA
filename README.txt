===============================================================================================
ParallaxBA: Bundle Adjustment with Parallax Angle Feature Parametrization
Version: 1.0
===============================================================================================

Copyright (C) 2013 Liang Zhao, Shoudong Huang, Yanbiao Sun and Gamini Dissanayake
University of Technology, Sydney, Australia

C/C++ sourse code for ParallaxBA: Bundle Adjustment with Parallax Angle Feature Parametrization

Authors:  Liang Zhao         -- Liang.Zhao-1@uts.edu.au 
          Shoudong Huang     -- Shoudong.Huang@uts.edu.au
          Yanbiao Sun        -- syb51@pku.edu.cn
	  Gamini Dissanayake -- Gamini.Dissanayake@uts.edu.au

          Centre for Autonomous Systems
          Faculty of Engineering and Information Technology
          University of Technology, Sydney
          NSW 2007, Australia
-----------------------------------------------------------------------------------------------
License
-----------------------------------------------------------------------------------------------

ParallaxBA: by Liang Zhao, Shoudong Huang, Yanbiao Sun, Gamini Dissanayake is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

1) You can freely use and modify this code.

2) If you want to distribute code based on this one, it has to be done under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

If you use this code for academic work, please reference:

      Liang Zhao, Shoudong Huang, Lei Yan, and Gamini Dissanayake,
      Parallax angle parametrization for monocular SLAM,
      IEEE International Conference on Robotics and Automation (ICRA), 2011.

      Liang Zhao, Shoudong Huang, Yanbiao Sun, Lei Yan and Gamini Dissanayake,
      ParallaxBA: Bundle Adjustment
using Parallax Angle Feature Parametrization,
      Submitted to International Journal of Robotics Research (IJRR), August 2013. 
     

3) For commercial uses, please contact the authors.

-----------------------------------------------------------------------------------------------
Quick start
-----------------------------------------------------------------------------------------------

On Windows platform, directly open the "demo.sln" solution and uncomment some codes corresponding
to one dataset in _tmain() function.

On Linux platform, install and run ParallaxBA as fowllows:

- sudo apt-get install cmake
- sudo apt-get install libeigen3-dev
- sudo apt-get install libsuitesparse-dev

- cd linux
- mkdir build
- cd build
- cmake ..
- make
- sudo make install

- cd Village
- ParallaxBA -cam Cam90.txt -fea Feature90.txt -calib cal90.txt -solve GN -report report.txt

-----------------------------------------------------------------------------------------------
Support
-----------------------------------------------------------------------------------------------

Questions, comments, suggestions, discussions and bug reports are welcomed. 

Please email to syb51@pku.edu.cn or Liang.Zhao-1@uts.edu.au.

