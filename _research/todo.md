---
title: To Do List
author: Aaron Kaw
---

# Stage 1: TL Field
1. Change medium range and depth to be a property of the source.
2. Investigate speed-up for field calculation:
   1. Mutating beam loop?
   2. Make initializing `Beam` method.
   3. Make mutating `Beam` method.
3. Make mature plotting function. Dispatch on:
   1. Single ray
   2. Vector of rays
   3. Single beam
   4. Vector of beams
   5. Pressure field (TL)
4. Implement more robust Gaussian beam width theory.
5. Repair matrix/vector inputs for range/depth dependent ocean parameters.
6. Investigate speed-up for ray trace (already pretty fast):
   1. Why is it still allocating?
7. General speed up:
   1. Use multithreading?
   2. Use GPU?
8. Implement progress bar.
9.  Test on a variety of scenarios:
   1. Clean up scenario script.
   2. Make plotting loops for fields.

# Stage 2: Sonar Equations
1. Detection index
2. Detection threshold
3. Probability of detection
4. Signal excess
   1. How to deal with zero detection index
5. 