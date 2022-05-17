#ifndef __TYPES_H__
#define __TYPES_H__

namespace NOL
{
#ifndef DATA_TYPE
#define DATA_TYPE
using FLOAT = double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
#endif

enum class LineSearchType
{
    GOLDENSECTION,
    ARMIJO,
    GOLDSTEIN,
    WOLFE,
    STRONGWOLFE
};

enum class QuasiNewtonSearchType
{
    SR1,
    DFP,
    BFGS
};

enum class ConjugateGradientType
{
    FR,
    PRP,
    PRP_PLUS,
    CD,
    DY
};

enum class TerminationCriterionType
{
    GK_NORM,
    DELTA_XK,
    DELTA_F
};

enum class OptimizationMethodType
{
    SD,
    NEWTON,
    DAMPED_NEWTON,
    LM,
    QUASI_NEWTON,
    ConjugateGradient
};

} // namespace NOL
#endif
