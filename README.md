# tmpl8rt_rust

## Building the Project
Ensure that you have `cargo` installed on your system. A guide to install `cargo` can be found [here](https://doc.rust-lang.org/cargo/getting-started/installation.html)

To build the project in Release mode run `cargo build --release`

To run the project in Release mode run `cargo run --release`

Running in debug mode is not recommended, but can be done by removing the `--release` flag for `cargo`

### Additional Linux Packages
The following packages need to be installed as well when you are running on Linux, as GLFW will not compile without it
- xorg-dev

## Controls
To control the project you can use `WASD` to move around and the arrow keys to look around

You can use the numbers `1`, `2`, `3` to get into the different viewpoints used in the experiment

You can use `C` to get the current result of the max intersection tests, minimum intersection tests, and average intersection tests. Furthermore, a screenshot is made and exported to `./output/`.

You can use `Space` to pause the rendering, which is useful for when the application gets stuck and you cannot use the GUI.

Other controls for ray tracing can be found in the GUI of imgui.
