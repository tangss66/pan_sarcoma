if (interactive() && Sys.getenv(TERM_PROGRAM) == vscode) {
  if (httpgd %in% .packages(all.available = TRUE)) {
    options(vsc.plot = FALSE)
    options(device = function(...) {
      httpgdhgd(silent = TRUE)
      .vsc.browser(httpgdhgd_url(history = FALSE), viewer = Beside)
    })
  }
}