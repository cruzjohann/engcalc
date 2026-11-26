# Transfer matrix checks against Chaudhry

Frequency response of pipe systems by the transfer matrix method, checked against
M. H. Chaudhry, *Applied Hydraulic Transients*, 3rd ed., Springer, 2014.

| File | What it does |
| --- | --- |
| `chaudhry_example_8_2.m` | Example 8-2 (reservoir, pipe, dead-end branch, pipe, oscillating valve). Runs in MATLAB and GNU Octave. |
| `chaudhry_example_8_2.py` | The same calculation in Python, plus the frequency response plot. |
| `appendix_d_check.py` | The two-pipe series system from Appendix D, compared row by row with the book's printed output. |
| `free_gas_example_8_2.py` | Example 8-2 with free gas in the line. Wave speed from Eq. 9-18 (Section 9-5), then the same frequency response. |

## Results

Example 8-2 at ωr = 2.0:

| Quantity | Book (hand calc) | Code |
| --- | --- | --- |
| hr = 2\|h\|/Ho | 0.170 | 0.1687 |
| qr = 2\|q\|/Qo | 0.390 | 0.3910 |
| Phase, head vs gate | −102.29° | −102.1765° |
| Phase, flow vs gate | −12.29° | −12.1765° |

The full sweep reproduces Fig. 8-17c. All 20 rows of the Appendix D output match to the
three decimals printed; the largest phase difference is 0.005°.

## Running

```
octave --eval "chaudhry_example_8_2"
python3 chaudhry_example_8_2.py
python3 appendix_d_check.py
python3 free_gas_example_8_2.py
```

The Python scripts need `numpy`; the plot needs `matplotlib`.

## Free gas

With 0.01% gas at atmospheric pressure, the wave speed in a rigid water pipe drops from 1481 to
833 m/s, which agrees with the book's statement in Section 2-6 that 1 part of air in 10,000 halves it.
In Example 8-2 at 10.8 bar abs, 0.1% gas moves the first resonance from 0.330 to 0.247 Hz, and the
head response at the book's hand-calculation frequency (0.667 Hz) goes from 0.169 to 0.492.

## Notes

- Example 8-2 uses g = 9.80665 m/s². The Appendix D program uses g = 9.81 m/s².
- Example 8-2 normalises head as 2|h|/Ho; Appendix D prints |h|/Ho.
- The wave speeds in both cases are the book's given values, not calculated.
- The method is linear and these cases are frictionless. It finds resonances, mode shapes and
  amplification for small oscillations. It does not give design pressures for large transients
  such as valve slam or pump trip, and it does not model cavitation.

The book's text and figures are not reproduced here; equations and figures are cited by number.

## Licence

MIT, see `LICENSE` in this folder.
