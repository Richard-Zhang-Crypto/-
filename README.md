# Cryptographic Mathematics Challenge - Problem 3 Source Code Repository

This repository contains the implementation code for Problem 3 of the Cryptographic Mathematics Challenge, covering the calculation of matrix order, linear hull search, and correlation estimation.

## 📁 Directory Structure & Function Descriptions

### 1. Matrix Order Calculation (Corresponding to Chapter 2 of the Paper)
* **`Calculate the order of the matrix.py`**
  * **Description**: Calculates the order of the matrix involved in the 32-bit SPN structure algorithm.

### 2. Linear Hull Search & Correlation Estimation (Corresponding to Chapter 3 of the Paper)
* **`linear_hull_search.cpp`**
  * **Description**: Searches for high-correlation linear hulls for a user-specified number of rounds.
* **`cor_approximate_search.cpp`**
  * **Description**: Estimates the correlation of the linear hulls provided in `test_data.txt`.
  * **Output Metrics**: Estimated values, true values, execution time, scores, and the number of successful approximations.

### 3. 8-Round Linear Hull Construction (Corresponding to Chapter 4 of the Paper)
* **`linear_hull_of_round.m`** (MATLAB)
  * **Description**: Outputs the 8-round linear hull constructed in Chapter 4 and its corresponding correlation.
  * **Output File**: Results are automatically saved to `result8.txt`.

---
