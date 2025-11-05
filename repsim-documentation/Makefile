# Makefile for Quarto documentation (run from docs/)
# Shortcuts for building and previewing

SHELL := /bin/bash
SITE_DIR := ../docs
QUARTO ?= /Applications/quarto/bin/quarto

.PHONY: all build preview clean rebuild open build_preview rebuild_preview check help

# --------------------------------------------------
# Core targets
# --------------------------------------------------

all build:
	@echo "==> Building site with $(QUARTO)..."
	@$(QUARTO) render

preview:
	@echo "==> Starting live preview with $(QUARTO)..."
	@$(QUARTO) preview

clean:
	@echo "==> Cleaning Quarto build artifacts..."
	@rm -rf "$(SITE_DIR)"
	@find . -name ".quarto" -type d -prune -exec rm -rf {} +

rebuild: clean build

open:
	@if [ -f "$(SITE_DIR)/index.html" ]; then \
		if command -v open >/dev/null 2>&1; then open "$(SITE_DIR)/index.html"; \
		elif command -v xdg-open >/dev/null 2>&1; then xdg-open "$(SITE_DIR)/index.html"; \
		else echo "Open $(SITE_DIR)/index.html in your browser."; fi; \
	else \
		echo "No site found. Run 'make' first."; exit 1; \
	fi

# --------------------------------------------------
# Combined workflow shortcuts
# --------------------------------------------------

build_preview:
	@echo "==> Building then launching preview..."
	@$(QUARTO) render && $(QUARTO) preview

rebuild_preview:
	@echo "==> Cleaning, rebuilding, then launching preview..."
	@rm -rf "$(SITE_DIR)" && \
		find . -name ".quarto" -type d -prune -exec rm -rf {} + && \
		$(QUARTO) render && $(QUARTO) preview

# --------------------------------------------------
# Diagnostics
# --------------------------------------------------

check:
	@echo "==> Quarto path : $(QUARTO)"
	@echo "==> Quarto version"
	@$(QUARTO) --version || true
	@echo "==> Python seen by reticulate (from R)"; \
	R -q -e 'suppressPackageStartupMessages({library(reticulate)}); try(py_discover_config(), silent=TRUE)'

help:
	@echo "Targets:"
	@echo "  make build             Build site (alias: make)"
	@echo "  make preview           Live preview"
	@echo "  make clean             Remove build artifacts"
	@echo "  make rebuild           Clean + build"
	@echo "  make open              Open built index.html"
	@echo "  make build_preview     Build once, then preview"
	@echo "  make rebuild_preview   Clean, build, then preview"
	@echo "  make check             Show Quarto path/version and interpreter info"
