#ifndef RESOURCES_HPP
#define RESOURCES_HPP

#include "data/data.hpp"

namespace Anaquin
{
    Scripts PlotNorm();
    Scripts PlotGene();
    Scripts PlotInsert();
    Scripts PlotLinear();
    Scripts PlotAllele();
    Scripts PlotSomatic();
    Scripts PlotLDensity();
    Scripts PlotKSomatic();
    Scripts PlotKSynthetic();
    Scripts PlotStrelkaROC();
    Scripts PlotQualFilter();
    Scripts PlotSplitSomatic();

    Scripts rna();
    Scripts meta();
    Scripts norm();
    Scripts split();
    Scripts somatic();
    Scripts broadBAM();
    Scripts broadVCF();
    Scripts germline();
    Scripts partition();
    Scripts calibrate();
    
    Scripts CSS();
    Scripts PReport();
    Scripts NormHTML();
    Scripts GermHTML();
    Scripts SomaHTML();
    Scripts GSplitHTML();
    Scripts RSplitHTML();
    Scripts MSplitHTML();
    Scripts CalibrateHTML();
}

#endif
