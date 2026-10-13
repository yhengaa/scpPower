# One-page export of the existing plausibility note; no model is fitted here.
note <- readLines("reports/M1/plausibility-note.md", warn = FALSE)
note <- gsub("\\[([^]]+)\\]\\([^)]*\\)", "\\1", note)
note <- gsub("[`*#]", "", note)
note <- unlist(lapply(note, function(line)
    if (!nzchar(line)) "" else strwrap(line, width = 114L)))
grDevices::pdf("reports/M1/plausibility-note.pdf", width = 8.27, height = 11.69,
               title = "scpPower M1 plausibility note: stopped calibration")
graphics::par(mar = c(0, 0, 0, 0))
graphics::plot.new(); graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1),
                                          xaxs = "i", yaxs = "i")
line_height <- min(0.0145, 0.90 / length(note))
for (i in seq_along(note)) graphics::text(0.06, 0.95-(i-1)*line_height,
    note[[i]], adj = c(0, 1), family = "Helvetica", cex = if (i == 1L) 1.0 else 0.66)
grDevices::dev.off()
