#include "sequins/genomics/g_calibrate.hpp"

using namespace Anaquin;

typedef GCalibrate::Options Options;

GCalibrate::Stats GCalibrate::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    GBroadBam::Options o_(cloneO(o)); o_.index = o.index; o_.debug = true; o_.showGen = false;
    o_.meth  = o.meth;
    o_.origW = o.work;
 
    if (o.debug)
    {
        switch (o_.meth)
        {
            case CalibrateMethod::None:         { o.logInfo("None");         break; }
            case CalibrateMethod::Mean:         { o.logInfo("Mean");         break; }
            case CalibrateMethod::Median:       { o.logInfo("Median");       break; }
            case CalibrateMethod::Percent:      { o.logInfo("Percent");      break; }
            case CalibrateMethod::SampleMean:   { o.logInfo("SampleMean");   break; }
            case CalibrateMethod::SampleMedian: { o.logInfo("SampleMedian"); break; }
        }
    }
    
    // Generate outputs like BroadBAM (we're going to move files later)
    GBroadBam::report(f1, f2, o_);
    
    auto m1 = [&](const FileName &src, const FileName &dst)
    {
        mv(o_.work + "/" + src, o.work + "/" + dst);
    };

    auto m2 = [&](const FileName &src, const FileName &dst)
    {
        mv(o.work + "/" + src, o.work + "/" + dst);
    };

    if (o.writeS) { m1("sample.bam", "sample.bam"); }
    if (o.writeD) { m1("trimmed.bam", "sequin.bam"); } // Not from the untrimmed "sequin.bam"
    if (o.writeC) { m1("sequin_calibrated.bam", "calibrated.bam"); }

    m1("merged.bam",    "merged.bam");
    m1("broad_bam.txt", "calibrate_report.txt");

    if (o.debug)
    {
        m1("broadBAM_regions.tsv",   "calibrate_regions.tsv");
        m1("broadBAM_variants.tsv",  "calibrate_variants.tsv");
        m1("broadBAM_synthetic.tsv", "calibrate_synthetic.tsv");
        m1("broadBAM_errors.tsv",    "calibrate_errors.tsv");
        m1("broadBAM_features.tsv",  "calibrate_features.tsv");
    }

    return GCalibrate::Stats();
}

static void writeReport(const SOptions &o)
{
    auto writeSomatic = [&](const FileName &file, const FileName &src)
    {
        o.generate(file);
        o.writer->open(file);
        o.writer->write((boost::format(PlotSomatic()) % date()
                                                      % src
                                                      % o.work
                                                      % src
                                                      % "Allele Frequency Ladder"
                                                      % "Expected Allele Frequency (log2)"
                                                      % "Measured Allele Frequency (log2)"
                                                      % "data$EXP_FREQ"
                                                      % "data$OBS_FREQ_CALIB").str());
        o.writer->close();
    };
}

void GCalibrate::report(const FileName &file, const Options &o)
{
    analyze(file, "", o);
    writeReport(o);
}

void GCalibrate::report(const FileName &f1, const FileName &f2, const Options &o)
{
    analyze(f1, f2, o);
    writeReport(o);
}
