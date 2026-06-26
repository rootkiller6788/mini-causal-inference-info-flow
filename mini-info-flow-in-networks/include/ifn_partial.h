#ifndef IFN_PARTIAL_H
#define IFN_PARTIAL_H
#include "ifn_core.h"

/* Partial Information Decomposition (PID).
 * Decomposes I(X;Y,Z) into:
 *   - Unique: info only from Y or only from Z
 *   - Redundant: info present in both Y and Z
 *   - Synergistic: info only available from Y and Z jointly */

typedef struct {
    double unique_y;
    double unique_z;
    double redundant;
    double synergistic;
    double total_mi;
} IFN_PIDResult;

IFN_PIDResult ifn_pid_decompose(const IFN_TimeSeries* target,
    const IFN_TimeSeries* source1, const IFN_TimeSeries* source2,
    int n_bins);

IFN_PIDResult ifn_pid_decompose_discrete(const int* target,
    const int* src1, const int* src2, int T, int n_states_t,
    int n_states_s1, int n_states_s2);

/* Specific PID measures */
double ifn_unique_information(const IFN_TimeSeries* target,
    const IFN_TimeSeries* source, const IFN_TimeSeries* other, int n_bins);

double ifn_redundant_information(const IFN_TimeSeries* target,
    const IFN_TimeSeries* src1, const IFN_TimeSeries* src2, int n_bins);

double ifn_synergistic_information(const IFN_TimeSeries* target,
    const IFN_TimeSeries* src1, const IFN_TimeSeries* src2, int n_bins);

/* Synergy from transfer entropy perspective */
double ifn_transfer_synergy(const IFN_TimeSeries* source1,
    const IFN_TimeSeries* source2, const IFN_TimeSeries* target,
    int k, int l);

void ifn_pid_print(const IFN_PIDResult* pid);

#endif
