# Deploy app.R to ShinyApps

rsconnect::deployApp(
  appName = "protAssay",
  appDir = here::here(),
  appFileManifest = here::here("dev", "app_file_manifest.txt")
)
