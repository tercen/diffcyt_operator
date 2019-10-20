# asinh operator

#### Description
`diffcyt` operator a differential analysis of flow cyto data and indicates which marker and cluster combinations are relevant.

##### Usage
Input projection|.
---|---
`col`    | `group_id`, `patient_id`
`row`    | `marker_class`, `marker_name`
`y-axis` | values representing measurement


Input parameters|.
---|---
`analysis_type` | can be either `DA` (Differential Abundance) or `DS` (Differential State)

Output relations|.
---|---
`cluster_id`| character, cluster  name, per cluster, DA and DS output
`LogFC`     | numeric, log fold change, per cluster, DA
`LR`        | numeric, lr, per cluster, DA output
`p_val`     | numeric, p value, per cluster, DA output
`p_adj`     | numeric, adjusted p value, per cluster, DA output
`LogFC`     | numeric, log fold change, per cluster-marker, DS output
`p_val`     | numeric, p value, per cluster-marker, DS output
`p_adj`     | numeric, adjusted p value, per cluster-marker, DS output
`t`         | numeric, t, per cluster-marker, DS output
`B`         | numeric, B, per cluster-marker, DS output
`AvgExp`    | numeric, adjusted p value, per cluster-marker, DS output

##### Details

Performs differential analysis (abundance or state). See the `diffcyt::diffcyt` function in the Bioconductor R pacakge.

#### References
see the github for documentation,
https://github.com/lmweber/diffcyt


##### See Also
[t-test](https://github.com/lmweber/ttest_operator), [anova](https://github.com/tercen/anova_operator), [rfImp](https://github.com/tercen/rfImp_operator)

#### Examples
