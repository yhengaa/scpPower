# Batch adjacent fixed points in one allocation to reduce scheduler/startup cost.
# Existing complete points are retained; each numerical seed stays unchanged.
args<-commandArgs(trailingOnly=TRUE);chunk<-as.integer(args[1])
rows<-seq.int((chunk-1L)*35L+1L,min(chunk*35L,945L))
for(i in rows) {
 path<-file.path('reports/M3/points',sprintf('%04d-manifest.json',i))
 if(file.exists(path)){cat('Retain completed grid point',i,'\n');next}
 Sys.setenv(SCPPOWER_GRID_TASK=as.character(i))
 source('inst/scripts/run_m3_grid.R',local=new.env(parent=globalenv()))
 invisible(gc())
}
