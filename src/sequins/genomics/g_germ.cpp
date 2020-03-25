#include <cmath>
#include "data/resources.hpp"
#include "writers/vcf_writer.hpp"
#include "writers/file_writer.hpp"
#include "sequins/genomics/g_germ.hpp"

using namespace Anaquin;

typedef GGerm::Stats Stats;
typedef GVariant::Options Options;

template <typename T, typename O> void writeFP(const std::vector<T> &x,
                                               const std::string &label,
                                               const std::string &format,
                                               const O &o)
{
    const auto &r = Standard::instance().gen;
    
    for (const auto &i : x)
    {
        if (i.rID.empty())
        {
            continue;
        }
        
        const auto name = (i.var && i.alt && i.ref ? i.var->name : MISSING);
        o.writer->write((boost::format(format) % i.rID
                                               % i.qry.cID
                                               % i.qry.l.start
                                               % label
                                               % gt2str(i.qry.gt)
                                               % var2str(i.qry.type())
                                               % (name != MISSING ? std::to_string(r.af(name)) : MISSING)
                                               % i.qry.dp(0)
                                               % i.qry.dp(1)
                                               % i.qry.obsAF()
                                               % S2(i.qry.qual[0])
                                               % MISSING
                                               % r.a1()->strForLocus(i.qry.cID, i.qry.l)).str());
    }
}

template <typename T, typename O> void writeGD(const std::string &format,
                                               const T &x,
                                               const O &o)
{
    writeFP(x.ds.fps, "FP", format, o);
}

template <typename T, typename O> void writeSV(const std::string &format,
                                               const T &x,
                                               const O &o)
{
    const auto &r = Standard::instance().gen;
    
    for (const auto &i : x.hs.vs)
    {
        o.writer->write((boost::format(format) % i.name
                                               % i.cID
                                               % i.l.start
                                               % "SV"
                                               % gt2str(i.gt)
                                               % var2str(i.type())
                                               % MISSING
                                               % i.dp(0)
                                               % i.dp(1)
                                               % i.obsAF()
                                               % S2(i.qual[0])
                                               % MISSING
                                               % r.a1()->strForLocus(i.cID, i.l)).str());
    }
}

template <typename T, typename O> void writeGQ(const std::string &format,
                                               const T &stats,
                                               const O &o)
{
    const auto &r = Standard::instance().gen;
    const auto r4 = r.r4()->ginters();
    const auto r5 = r.r5()->ginters();

    for (const auto &i : r.v1()->data.vars())
    {
        if (isGermA(i.name))
        {
            // Can we find this sequin?
            const auto isTP = stats.ds.findTP(i.name);
            
            // Locus for the sequin
            const auto l5 = r5.at(i.cID).overlap(i.l)->l();
            
            if (isTP)
            {
                // Called variant (if found)
                const auto &q = isTP->qry;

                o.writer->write((boost::format(format) % i.name
                                                       % i.cID
                                                       % i.l.start
                                                       % "TP"
                                                       % gt2str(i.gt)
                                                       % var2str(i.type())
                                                       % r.af(i.name)
                                                       % q.dp(0)
                                                       % q.dp(1)
                                                       % q.obsAF()
                                                       % S2(q.qual[0])
                                                       % r4.at(i.cID).length(l5)
                                                       % r.a1()->strForLocus(i.cID, i.l)).str());
            }

            // Failed to detect the variant
            else
            {
                o.writer->write((boost::format(format) % i.name
                                                       % i.cID
                                                       % i.l.start
                                                       % "FN"
                                                       % gt2str(i.gt)
                                                       % var2str(i.type())
                                                       % r.af(i.name)
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % MISSING
                                                       % r4.at(i.cID).length(l5)
                                                       % r.a1()->strForLocus(i.cID, i.l)).str());
            }
        }
    }
}

template <typename T, typename O> void writeQuin(const FileName &file, const T &stats, const O &o)
{
    const auto &r = Standard::instance().gen;
    const auto format = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%\t%11%\t%12%%13%";
    
    o.writer->open(file);
    o.writer->write((boost::format(format) % "NAME"
                                           % "CHROM"
                                           % "POSITION"
                                           % "LABEL"
                                           % "GENOTYPE"
                                           % "TYPE"
                                           % "EXP_FREQ"
                                           % "REF_DEPTH"
                                           % "VAR_DEPTH"
                                           % "OBS_FREQ"
                                           % "QUAL"
                                           % "SIZE"
                                           % r.a1()->strForKeys()).str());

    writeGQ(format, stats, o); // Write true positives and false negatives
    writeGD(format, stats, o); // Write false positives
    writeSV(format, stats, o); // Write sample variants
    
    o.writer->close();
}

static void writeTable(const FileName &file, const GGerm::Stats &stats, const GVariant::Options &o)
{
    const auto f = "%1%\t%2%\t%3%\t%4%\t%5%\t%6%\t%7%\t%8%\t%9%\t%10%";
    
    o.generate(file);
    o.writer->open(file);
    o.writer->write("VARIANTS\tTOTAL\tTP\tFP\tFN\tSN\tPC\tFP_KB\tDEPTH\tSAMPLE");

    auto write = [&](const Label &name, const GVariant::TableRow &x)
    {
        o.writer->write((boost::format(f) % name
                                          % (x.valid ? replaceNA(x.nr)     : MISSING)
                                          % (x.valid ? replaceNA(x.tp)     : MISSING)
                                          % (x.valid ? replaceNA(x.fp)     : MISSING)
                                          % (x.valid ? replaceNA(x.fn)     : MISSING)
                                          % (x.valid ? replaceNA(x.sn())   : MISSING)
                                          % (x.valid ? replaceNA(x.pc())   : MISSING)
                                          % (x.valid ? replaceNA(x.fpKB()) : MISSING)
                                          % (x.valid ? replaceNA(x.depth)  : MISSING)
                                          % (x.valid ? replaceNA(x.sample) : MISSING)).str());
    };

    write("Total_variants",              GVariant::getTotal(o));
    write("SNV_type",                    GVariant::getSNV(o));
    write("Indel_type",                  GVariant::getIndel(o));
    write("Homozygous_Genotype",         GVariant::getHom(o));
    write("Heterozygous_Genotype",       GVariant::getHet(o));
    write("CodingRegion_GeneContext",    GVariant::getCode(o));
    write("NoncodingRegion_GeneContext", GVariant::getNCode(o));
    // TODO: write("ATrich_GCcontent",            GVariant::getAT(o));
    // TODO: write("GCrich_GCcontent",            GVariant::getGC(o));
    write("SINE_MobileElement",          GVariant::getSine(o));
    write("LINE_MobileElement",          GVariant::getLine(o));
    write("LTR_MobileElement",           GVariant::getLTR(o));
    write("DNA_repeat_MobileElement",    GVariant::getDNA(o));
    write("Mono_SimpleRepeat",           GVariant::getDi(o));
    write("Di_SimpleRepeat",             GVariant::getTri(o));
    write("Tri_SimpleRepeat",            GVariant::getMono(o));
    write("Quad_SimpleRepeat",           GVariant::getQuad(o));

    o.writer->close();
}

static void writeSummary(const FileName &file,
                         const FileName &f1,
                         const FileName &f2,
                         const Stats &stats,
                         const GVariant::Options &o)
{
    const auto &r = Standard::instance().gen;

    // Total sequin size for this analysis
    const auto sl = GVariant::countSizeForG(o);
    
    const auto f = "1. ANALYSIS\n"
                   "Date:                            %1%\n"
                   "Anaquin version:                 %2%\n"
                   "Command:                         %3%\n"
                   "Reference variants:              %4%\n"
                   "Reference regions:               %5%\n"
                   "Reference stratification:        %6%\n"
                   "Input variant file:              %7%\n"
                   "TP sequin variants:              %8%\n"
                   "FP sequin variants:              %9%\n"
                   "FN sequin variants:              %10%\n\n"
                   "2. GERMLINE VARIANTS\n"
                   "Detected sequin variants (TP):   %11%\n"
                   "Undetected sequin variants (FN): %12%\n"
                   "Erroneous sequin variants (FP):  %13%\n"
                   "Sensitivity (SN):                %14%\n"
                   "Precision (PC):                  %15%\n"
                   "TP median depth:                 %16%\n"
                   "TP median QUAL score (TPMQ):     %17%\n"
                   "FP median depth:                 %18%\n"
                   "FP median QUAL score (FPMQ):     %19%\n"
                   "Analysed region size:            %20%\n"
                   "False positives per KB:          %21%\n\n"
                   "3. SNVs                          TP; FN; FP; SN; PC; TPMQ; FPMQ\n"
                   "Total:                           %22%\n"
                   "Homozygous:                      %23%\n"
                   "Heterozygous:                    %24%\n"
                   "High GC (>70%%):                  %25%\n"
                   "Moderate GC (30-70%%):            %26%\n"
                   "Low GC (<30%%):                   %27%\n"
                   "High confidence:                 %28%\n"
                   "Low confidence:                  %29%\n"
                   "Coding:                          %30%\n"
                   "Noncoding:                       %31%\n"
                   "ACMG genes:                      %32%\n"
                   "Pharmacogenes:                   %33%\n"
                   "Simple repeats:                  %34%\n"
                   "5-12nt homopolymer:              %35%\n"
                   ">12nt homopolymer:               %36%\n\n"
                   "4. INDELS                        TP; FN; FP; SN; PC; TPMQ; FPMQ\n"
                   "Total:                           %37%\n"
                   "Homozygous:                      %38%\n"
                   "Heterozygous:                    %39%\n"
                   "High GC (>70%%):                  %40%\n"
                   "Moderate GC (30-70%%):            %41%\n"
                   "Low GC (<30%%):                   %42%\n"
                   "High confidence:                 %43%\n"
                   "Low confidence:                  %44%\n"
                   "Coding:                          %45%\n"
                   "Noncoding:                       %46%\n"
                   "ACMG genes:                      %47%\n"
                   "Pharmacogenes:                   %48%\n"
                   "Simple repeats:                  %49%\n"
                   "5-12nt homopolymer:              %50%\n"
                   ">12nt homopolymer:               %51%";
    
    auto S = [&](const GVariant::TableRow &x)
    {
        return !x.valid ? MISSING : x.broad();
    };

    const auto tpq = GVariant::getTP(o).qs;
    const auto tpd = GVariant::getTP(o).ds;
    const auto fpq = GVariant::getFP(o).qs;
    const auto fpd = GVariant::getFP(o).ds;

    o.generate(file);
    o.writer->open(file);
    o.writer->write((boost::format(f) % date()
                                      % o.version
                                      % o.cmd
                                      % r.v1()->src
                                      % r.r3()->src
                                      % r.a1()->src
                                      % (f1 + (f1.empty() ? "" : " and " + f2))
                                      % (o.work + "/sequins_TP.vcf")
                                      % (o.work + "/sequins_FP.vcf")
                                      % (o.work + "/sequins_FN.vcf")
                                      % stats.ds.oc.tp()     // 11
                                      % stats.ds.oc.fn()     // 12
                                      % stats.ds.oc.fp()     // 13
                                      % S3(stats.ds.oc.sn()) // 14
                                      % S3(stats.ds.oc.pc()) // 15
                                      % (tpd.empty() ? MISSING : S2(med(tpd))) // 16
                                      % (tpq.empty() ? MISSING : S2(med(tpq))) // 17
                                      % (fpd.empty() ? MISSING : S2(med(fpd))) // 18
                                      % (fpq.empty() ? MISSING : S2(med(fpq))) // 19
                                      % sl                                     // 20
                                      % S4(1000.0 * ((Proportion) stats.ds.oc.fp() / sl)) // 21
                                      % S(GVariant::getSNV(o))                  // 22
                                      % S(GVariant::getHom(o, "SNP"))           // 23
                                      % S(GVariant::getHet(o, "SNP"))           // 24
                                      % S(GVariant::getHighGC(o, "SNP"))        // 25
                                      % S(GVariant::getModGC(o, "SNP"))         // 26
                                      % S(GVariant::getLowGC(o, "SNP"))         // 27
                                      % S(GVariant::getHighConf(o, "SNP"))      // 28
                                      % S(GVariant::getLowConf(o, "SNP"))       // 29
                                      % S(GVariant::getCode(o, "SNP"))          // 30
                                      % S(GVariant::getNCode(o, "SNP"))         // 31
                                      % S(GVariant::getACMG(o, "SNP"))          // 32
                                      % S(GVariant::getPGX(o, "SNP"))           // 33
                                      % S(GVariant::getRepeat(o, "SNP"))        // 34
                                      % S(GVariant::getShortRepeat(o, "SNP"))   // 35
                                      % S(GVariant::getLongRepeat(o, "SNP"))    // 36
                                      % S(GVariant::getIndel(o))                // 37
                                      % S(GVariant::getHom(o, "Indel"))         // 38
                                      % S(GVariant::getHet(o, "Indel"))         // 39
                                      % S(GVariant::getHighGC(o, "Indel"))      // 40
                                      % S(GVariant::getModGC(o, "Indel"))       // 41
                                      % S(GVariant::getLowGC(o, "Indel"))       // 42
                                      % S(GVariant::getHighConf(o, "Indel"))    // 43
                                      % S(GVariant::getLowConf(o, "Indel"))     // 44
                                      % S(GVariant::getCode(o, "Indel"))        // 45
                                      % S(GVariant::getNCode(o, "Indel"))       // 46
                                      % S(GVariant::getACMG(o, "Indel"))        // 47
                                      % S(GVariant::getPGX(o, "Indel"))         // 48
                                      % S(GVariant::getRepeat(o, "Indel"))      // 49
                                      % S(GVariant::getShortRepeat(o, "Indel")) // 50
                                      % S(GVariant::getLongRepeat(o, "Indel"))  // 51
                     ).str());
    o.writer->close();
}

GGerm::Stats GGerm::analyze(const FileName &f1, const FileName &f2, const Options &o)
{
    const auto &r = Standard::instance().gen;

    auto o_ = o;
    o_.es = tmpFile();
    o_.tp = o.work + "/sequins_TP.vcf";
    o_.fp = o.work + "/sequins_FP.vcf";
    o_.fn = o.work + "/sequins_FN.vcf";
    o_.tsvS = tmpFile();
    
    Stats stats;
    GVariant::analyze(stats, GGerm(), f1, f2, o_, Standard::instance().gen.v2());

    return stats;
}

void GGerm::report(const FileName &f1, const FileName &f2, const Options &o)
{
    const auto stats = analyze(f1, f2, o);

    auto o1 = o.debug ? o : cloneO(o);
    o1.tsvS = o1.work + "/germline_sequin.tsv";
    o1.writer = std::shared_ptr<Writer<>>(new FileWriter(o1.work));

    // Generating germline_sequin.tsv
    writeQuin("germline_sequin.tsv", stats, o1);

    if (o.debug)
    {
        // Generating germline_variant_table.tsv
        writeTable("germline_variant_table.tsv", stats, o1);
    }
    
    auto o2 = o;
    o2.tsvS = o1.tsvS;
    
    // Generating germline_report.txt
    writeSummary("germline_report.txt", f1, f2, stats, o2);
}
