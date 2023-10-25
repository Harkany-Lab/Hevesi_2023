npr <- c('Oxtr',  'Avpr1a', 'Avpr1b', 'Avpr2',
         'Tacr1', 'Mc1r', 'Mc3r', 'Mc4r', 'Oprl1',  'Tacr3',
         'Gpr83', 'Galr1', 'Npy1r', 'Npy2r', 'Npy4r', 'Npy4r2', 'Npy5r',
         'Sstr1', 'Sstr2', 'Sstr3', 'Mchr1',  'Oprd1', 'Oprk1', 'Oprm1',
         'Trhr', 'Hcrtr1', 'Hcrtr2', 'Qrfpr',   'Npffr1', 'Npffr2',
         'Prlhr',  'Ghr', 'Ghrhr',  'Grpr', 'Vipr1', 'Vipr2', 'Prokr2',
         'Nmur2',  'Nmur1', 'Nmbr',  'Kiss1r', 'Crhr1', 'Crhr2',
         'Cntfr', 'Cckar', 'Cckbr', 'Galr3')

np <- c('Adcyap1', 'Oxt', 'Avp', 'Tac1', 'Pomc', 'Pnoc', 'Tac2',
        'Nts', 'Gal', 'Agrp', 'Npy', 'Sst', 'Cartpt', 'Pmch',
        'Reln', 'Rxfp1', 'Penk', 'Pdyn', 'Trh', 'Hcrt', 'Qrfp',
        'Npw', 'Npvf', 'Ghrh', 'Grp', 'Vip', 'Nms', 'Nmu', 'Nmb',
        'Kiss1', 'Crh', 'Bdnf', 'Cntf', 'Cck')

neurotrans <- c("Slc17a6", "Slc17a7", "Slc17a8", "Slc1a1", "Slc1a2", "Slc1a6",
                "Gad1", "Slc32a1", "Slc6a1")
glut       <- c("Slc17a6", "Slc17a7", "Slc17a8", "Slc1a1", "Slc1a2", "Slc1a6")
glutr      <- c("Gria1", "Gria2", "Gria3", "Gria4", # iGlu AMPA receptors
                "Grid1", "Grid2", # iGlu delta receptors
                "Grik1", "Grik2", "Grik3", "Grik4", "Grik5", # iGlu kainate receptors
                "Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d", "Grin3a", "Grin3b", # iGlu NMDA receptors
                "Grm1", "Grm5", # mGluRs 1
                "Grm2", "Grm3", # mGluRs 2
                "Grm4", "Grm6", "Grm7", "Grm8"# mGluRs 3
)
gaba       <- c("Gad1", "Gad2", "Slc32a1", "Slc6a1")
gabar <- c(
  "Gabra1", "Gabra2", "Gabra3", "Gabra4", "Gabra5", "Gabra6",
  "Gabrb1", "Gabrb2", "Gabrb3",
  "Gabrg1", "Gabrg2", "Gabrg3",
  "Gabrd", "Gabre", "Gabrp", "Gabrq",
  "Gabrr1", "Gabrr2", "Gabrr3",
  "Gabbr1", "Gabbr2"
)

dopam <- c("Th", "Slc6a3", "Slc18a2", "Ddc", "Slc18a3")
ach <- c("Chat", "Slc18a3", "Ache", "Slc5a7")
genes.zh   <- c("Rbfox3", "Gja1", "Aqp4", "Slc1a3", "Gfap", "Slc17a6", "Slc17a7",
                "Gad1", "Gad2", "Slc18a2", "Plcb1", "Prkaca", "Adcy1", "Grin1", "Gap43",
                "Galr1", "Galr2", "Galr3", "Lmx1b", "Bax", "Dcaf5", "Ipo11", "Ntn1",
                "Slit1", "Slit2", "Robo1", "Dcc", "Nrg1", "Lmo4", "Syngap1",
                "Nf1", "Grm5", "Rfx3", "Ache", "Sox2", "Sox4", "Sox5",
                "Lamp", "Cad6", "Cad8", "Cad11", "Efna5", "Epha3", "Epha8",
                "Ephb2", "Ephb3", "Maoa", "Htr1b", "Slc6a4",
                "Jam2", "Galnt14", "Hs6st3", "Zfhx4", "Cck", "Lef1", "Ebf1",
                "Pou2f2", "Vgf", "Vegfc", "Gal")

gene_int <- c(npr, np, neurotrans, genes.zh) %>% unique()

