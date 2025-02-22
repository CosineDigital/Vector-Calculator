# Vector Calculator
![Vector Calculator App](https://raw.githubusercontent.com/CosineDigital/Vector-Calculator/master/image1.png)
<center>Vector Calculator - a powerful and fast calculator for engineering students</center>

## Features
- Robust Copy-Paste system for easy and fast calculations
- Position Vector from P to Q and magnitude
- Cross Product and magnitude
- Dot Product
- Distance from a point P to a line L with point L0
- Unit Vector and magnitude
- Scaled Unit Vector
- Cartesian Cosine Angles
- Cartesian to Spherical Coordinates
- Cartesian to Cylindrical Coordinates
- All values rounded to 3 decimal places
- Angles given in degrees
- Easy-to-switch themes via [Dear ImGui](https://github.com/ocornut/imgui)

## Notes
- Source code for the Vector Calculator can be found in [WalnutApp/src/WalnutApp.cpp](https://github.com/CosineDigital/Vector-Calculator/blob/master/WalnutApp/src/WalnutApp.cpp)

## Acknowledgements
- Created using [TheCherno/Walnut](https://github.com/TheCherno/Walnut)
- More information below

# Walnut

Walnut is a simple application framework built with Dear ImGui and designed to be used with Vulkan - basically this means you can seemlessly blend real-time Vulkan rendering with a great UI library to build desktop applications. The plan is to expand Walnut to include common utilities to make immediate-mode desktop apps and simple Vulkan applications.

Currently supports Windows - with macOS and Linux support planned. Setup scripts support Visual Studio 2022 by default.

![WalnutExample](https://hazelengine.com/images/ForestLauncherScreenshot.jpg)
_<center>Forest Launcher - an application made with Walnut</center>_

## Requirements
- [Visual Studio 2022](https://visualstudio.com) (not strictly required, however included setup scripts only support this)
- [Vulkan SDK](https://vulkan.lunarg.com/sdk/home#windows) (preferably a recent version)

## Getting Started
Once you've cloned, run `scripts/Setup.bat` to generate Visual Studio 2022 solution/project files. Once you've opened the solution, you can run the WalnutApp project to see a basic example (code in `WalnutApp.cpp`). I recommend modifying that WalnutApp project to create your own application, as everything should be setup and ready to go.

### 3rd party libaries
- [Dear ImGui](https://github.com/ocornut/imgui)
- [GLFW](https://github.com/glfw/glfw)
- [stb_image](https://github.com/nothings/stb)
- [GLM](https://github.com/g-truc/glm) (included for convenience)

### Additional
- Walnut uses the [Roboto](https://fonts.google.com/specimen/Roboto) font ([Apache License, Version 2.0](https://www.apache.org/licenses/LICENSE-2.0))
