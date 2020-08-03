# Align trajectory

This codes can be used to align trajectory. And show trajectory, show position after alignment, etc.

- The trajectory can be automatically aligned according to the timestamp.
- The alignment uses [umeyama algorithm](https://pdfs.semanticscholar.org/d107/231cce2676dbeea87e00bb0c587c280b9c53.pdf?_ga=2.264495439.1181657306.1595240335-198766482.1595240335).
- Trajectory format: timestamp x y z q_x q_y q_z q_w

![](./data/gui.gif)

## Dependencies
- Pangolin
- EIgen

## Usage

Compile: 

```bash
cd align_trajectory
mkdir build
cd build
cmake ..
make 
```

run:

```bash
cd align_trajectory
./align_trajectory
```

