from ruffus import *
import CGATPipelines.Pipeline as P
import os
import sys
import CGATPipelines.PipelineGeneInfo as PipelineGeneInfo
import CGAT.IOTools as IOTools

PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

example_pw = PARAMS['my_gene_info_pathway'].split(",")[0]
if example_pw == "all":
    example_pw = 'kegg'

example_hg = [str(host) for host in
              str(PARAMS['my_gene_info_homologene']).split(",")][0]
if example_hg == "all":
    example_hg = '10090'

mgiannotations = PARAMS['my_gene_info_annotations'].split(",")


@originate('allgenes.tsv')
def GetAndTranslateAllGenes(outfile):
    GeneAnnot = PipelineGeneInfo.EntrezGeneAnnotation(
        PARAMS['db_name'], PARAMS['entrez_email'])
    genelist = GeneAnnot.download_all(PARAMS['entrez_host'])

    Ens = PipelineGeneInfo.EnsemblAnnotation(PARAMS['my_gene_info_source'],
                                             PARAMS['db_name'])
    PipelineGeneInfo.runall(Ens, genelist, ['ensembl'], submit=True)
#    Ens.runall(genelist, ['ensembl'])

    Sym = PipelineGeneInfo.SymbolAnnotation(PARAMS['my_gene_info_source'],
                                            PARAMS['db_name'],
                                            PARAMS['entrez_host'],
                                            PARAMS['entrez_sciname'])
    PipelineGeneInfo.runall(Sym, genelist, ['symbol'], submit=True)
#    Sym.runall(genelist, ['symbol'])
    outf = IOTools.openFile(outfile, "w")
    for gene in genelist:
        outf.write("%s\n" % gene)
    outf.close()


@active_if('go' in mgiannotations)
@transform(GetAndTranslateAllGenes, suffix(".tsv"),
           'ensemblg2go$annot.load')
def AnnotateWithGO(infile, outfile):
    genelist = PipelineGeneInfo.readGeneList(infile)
    GO = PipelineGeneInfo.GoAnnotation(PARAMS['my_gene_info_source'],
                                       PARAMS['db_name'],
                                       PARAMS['my_gene_info_go'])
    PipelineGeneInfo.runall(GO, genelist, ['go'], submit=True)
#    GO.runall(genelist, ['go'])
    ont = PipelineGeneInfo.OntologyAnnotation('go',
                                              PARAMS['my_gene_info_goont'],
                                              PARAMS['db_name'])
#    ont.runall(genelist)
    PipelineGeneInfo.runall(ont, genelist, submit=True)

@active_if('pathway' in mgiannotations)
@transform(GetAndTranslateAllGenes, suffix(".tsv"),
           'ensemblg2%s$annot.load' % example_pw)
def AnnotateWithPathway(infile, outfile):
    genelist = PipelineGeneInfo.readGeneList(infile)
    PW = PipelineGeneInfo.PathwayAnnotation(PARAMS['my_gene_info_source'],
                                            PARAMS['db_name'],
                                            PARAMS['my_gene_info_pathway'])
   # PW.runall(genelist, ['pathway'])
    PipelineGeneInfo.runall(PW, genelist, ['pathway'], submit=True)


@active_if('homologene' in mgiannotations)
@transform(GetAndTranslateAllGenes, suffix(".tsv"),
           'ensemblg2symbol_Homo_sapiens$annot.load')
def AnnotateWithHomologene(infile, outfile):
    genelist = PipelineGeneInfo.readGeneList(infile)
    HG = PipelineGeneInfo.HomologeneAnnotation(PARAMS['my_gene_info_source'],
                                               PARAMS['db_name'],
                                               PARAMS[
                                                   'my_gene_info_homologene'],
                                               PARAMS['entrez_host'],
                                               PARAMS['entrez_email'])
#    HG.runall(genelist, ['homologene'])
    PipelineGeneInfo.runall(HG, genelist, ['homologene'], submit=True)

@follows(AnnotateWithHomologene)
@active_if(int(PARAMS['homologues_mgi']) == 1)
@transform('ensemblg2symbol_Mus_musculus$geneid.load', suffix(".load"),
           'ensemblg2mousepathway$annot.load')
def AnnotateWithMousePathway(infile, outfile):
    genelist = PipelineGeneInfo.getSymbols(infile)
    MP = PipelineGeneInfo.MousePathwayAnnotation(
        PARAMS['homologues_mousemine'],
        PARAMS['db_name'])
#    MP.runall(genelist)
    PipelineGeneInfo.runall(MP, genelist, submit=True)

@follows(AnnotateWithHomologene)
@active_if(int(PARAMS['homologues_mousepathway']) == 1)
@transform('ensemblg2symbol_Mus_musculus$geneid.load', suffix(".load"),
           'ensemblg2mgi$annot.load')
def AnnotateWithMGI(infile, outfile):
    genelist = PipelineGeneInfo.getSymbols(infile)
    MGI = PipelineGeneInfo.MGIAnnotation(
        PARAMS['homologues_mousemine'],
        PARAMS['db_name'])
#    MGI.runall(genelist)
    PipelineGeneInfo.runall(MGI, genelist, submit=True)

@follows(AnnotateWithHomologene)
@active_if(int(PARAMS['homologues_hpo']) == 1)
@transform('ensemblg2symbol_Homo_sapiens$geneid.load', suffix(".load"),
           'ensemblg2hpo$annot.load')
def AnnotateWithHPO(infile, outfile):
    genelist = PipelineGeneInfo.getSymbols(infile)
    HPO = PipelineGeneInfo.HPOAnnotation(
        PARAMS['homologues_humanmine'],
        PARAMS['db_name'])
#    HPO.runall(genelist)
    PipelineGeneInfo.runall(HPO, genelist, submit=True)
    ont = PipelineGeneInfo.OntologyAnnotation('hpo',
                                              PARAMS['homologues_hpoont'],
                                              PARAMS['db_name'])
    ont.runall(genelist)


@follows(AnnotateWithGO)
@follows(AnnotateWithPathway)
@follows(AnnotateWithHomologene)
@follows(AnnotateWithMGI)
@follows(AnnotateWithMousePathway)
@follows(AnnotateWithHPO)
@transform("allgenes.tsv", suffix(".tsv"), '.tsv')
def full(infile, outfile):
    pass

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
