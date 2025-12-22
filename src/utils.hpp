
#include <TMB.hpp>

// GMRF for random deficient componts
// 
// Given
// x = (x_A, x_B)^T
// where
// x_A has V_A that is full rank, and some are fixed
// x_B has V_B that has no rank (V_B = 0), and none are fixed
//
// Define
// P = | P_A,  P_AB |
//     | P_BA, P_B  |
//
// V = | V_A,  V_AB |
//     | V_BA, V_B  |
//
// M = | I-P_A,  P_AB  |  =  | M_A,  M_AB |
//     | P_B_A,  I-P_B |     | M_BA, M_B  |
//
// Calculate
// C = M_AB M_B^-1
// so
// C^T = (M_B^T)^-1 M_AB^T
// and
// Mtilda_A = M_A - M_AB M_BB^-1 M_BA
// Vtilda_A = V_A + C V_B C^T + C V_BA + V_AB C^T 
//          = V_A + C V_BA + V_AB C^T   (because V_B = 0)
// Q_A = Mtilda_A^T V_A^-1 Mtilda_A
//
// Then:
// x_A ~ GMRF( Q_A )
// mu_B = -M_BB^-1 M_BA x_A (conditional krigging)
//
// And 
// x_B = mu_B
// Because 
// x_B ~ MVN( mu_B, Q_B^-1 )
// And:
// Q_B = M_BB^T V_B^-1 M_BB 
// so 
// Q_B^-1 = 0 (because V_B = 0)

