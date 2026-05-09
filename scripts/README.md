# Visualization

Visualizations are computed using Python with Matplotlib animation.

```bash
pip install numpy matplotlib
```

Then run the script with either a csv file (*t x1 v1 x2 v2 x3 v3* no header)

```bash
python3 scripts/visualize.py --file output.csv
```

or directly from `simulate`

```bash
./simulate X1 X2 X3 X4 X5 X6 X7 X8 | python3 scripts/visualize.py
```

Save the animation

```bash
python3 scripts/visualize.py --save output.gif
```