# Make org.db database for Xenopus tropicalis
# Run once
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("AnnotationHub", quietly = TRUE)){
  BiocManager::install("AnnotationHub")
}

if (!require("AnnotationForge", quietly = TRUE)){
  BiocManager::install("AnnotationForge")
}

hub <- AnnotationHub()
q <- query(hub, c("org.Xenopus_tropicalis.eg.sqlite","Xenopus tropicalis"))


##### Build OrgDb

file.copy(AnnotationHub::cache(q[1]), "org.xenTro.eg.sqlite")


seed <- new("AnnDbPkgSeed", 
            Package = "org.xenTro.eg.db", 
            Version = "0.0.1",
            Author = "Mauricio Roza", 
            Maintainer = "Mauricio Roza <mauricio.roza@aces.su.se>", 
            PkgTemplate = "NOSCHEMA.DB", 
            AnnObjPrefix = "org.xenTro.eg", 
            organism = "Xenopus tropicalis", 
            species = "Xenopus tropicalis", 
            biocViews = "annotation", 
            manufacturerUrl = "none", 
            manufacturer = "none", 
            chipName = "none")

makeAnnDbPkg(seed, dbfile = "org.xenTro.eg.sqlite")

install.packages("org.xenTro.eg.db/", type = "source", repos = NULL)

rm(hub)
rm(q)
rm(seed)
gc()

unlink("org.xenTro.eg.db/", recursive = TRUE)

unlink("org.xenTro.eg.sqlite")
