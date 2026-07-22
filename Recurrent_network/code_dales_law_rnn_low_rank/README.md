# Low-Rank Recurrent BANFF Experiments

This folder contains the recurrent MATLAB implementation associated with the
manuscript. The recurrent state is decoded through a fixed low-rank scaffold,
and only task-specific neuronal bias currents are trained.

## Model

For encoder matrix `eta`, decoder matrix `phi`, and bias vector `bias`, the
decoded state obeys

```matlab
dx/dt = -x + phi' * tanh(eta*x + bias)
```

At sampled states, training therefore supervises the static network output
against `x + dx/dt`. `bias_gradient_descent.m` changes only `bias`; `eta` and
`phi` remain fixed. Its momentum coefficient is `0.9`.

## Run Order

Set the MATLAB current folder to this directory before running these scripts.
The spreadsheet is resolved relative to the script, but generated MAT files and
figures are written to the current folder.

### Standard multisystem network: Figures 1I-K, 2, S2, and S3

Run:

```matlab
rnn_learning
```

`net_config.m` creates the shared fixed network:

- neurons: `N = 6000`
- decoded rank: `3`
- encoder values: `-4`, `0`, or `4` with probabilities `0.25`, `0.5`, `0.25`
- decoder distribution: Gaussian with standard deviation `1/sqrt(N)`
- random seed: `1` using the Mersenne Twister

`rnn_learning.m` uses `1000` state samples, `100` validation samples,
learning rate `0.1`, momentum `0.9`, and `300000` updates. It trains the first
40 spreadsheet systems by default, matching the retained experiment script.
Edit `SYSTEM_ROWS` only when intentionally selecting a subset.

### Multiple Lorenz bias solutions: Figure 3

Run:

```matlab
rnn_lorenz_bias_ensemble
```

The script reuses the standard 6000-neuron network and trains 20 Lorenz bias
vectors. Each starts from an independent Gaussian vector with standard
deviation `3` and is trained for `10000` updates with learning rate `0.1` and
momentum `0.9`. It then performs both manuscript switching protocols:

- the first five bias vectors are switched every `400` time units during a
  `2000`-time-unit `ode45` rollout;
- one of all 20 vectors is selected at every forward-Euler step with
  `dt = 0.01` during the rapid-switching rollout.

The saved output is `lorenz_bias_ensemble_results.mat`.

### Dale-constrained network: Supplementary Figure S4

Run:

```matlab
rnn_learning_dales_law
```

`net_config_dales_law.m` creates 3000 excitatory and 3000 inhibitory neurons.
Encoder entries are independently `0` or `4` with equal probability. Decoder
rows carry the sign of the presynaptic neuron and use rounded Gaussian
components scaled by `1/sqrt(N)`. Biases start from the standard normal
distribution. Training uses learning rate `0.02`, momentum `0.9`, and `30000`
updates. The six displayed spreadsheet rows are exposed as `SYSTEM_ROWS`.

## Source Files

- `bias_gradient_descent.m`: shared bias-only gradient descent with momentum.
- `net_int.m`: low-rank decoded-state ODE right-hand side.
- `net_int_large.m`: full recurrent-state ODE helper.
- `int_dyn.m`: target-system vector fields.
- `dynamics_list.xlsx`: system labels, dimensions, inputs, initial conditions,
  and scaling flags.
- `dynamics_list.m`: legacy textual task list retained for compatibility.
- `readme.docx`: original notes; this Markdown file is the current guide.

## Reproducibility and Runtime

The scripts set `rng(1,'twister')` before creating the fixed network and sample
sets. Re-running from a clean MATLAB session preserves the explicit random
sequence. Exact trajectories can still vary slightly across MATLAB versions or
platforms because `ode45` is adaptive. The standard network is large and the
300000-update multisystem run is computationally expensive; test a small
temporary value of `SYSTEM_ROWS` before launching a complete reproduction.
