#!/usr/bin/env python3
"""Static regression checks for the version-controlled Apptainer recipe."""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
RECIPE = ROOT / "container" / "sedimix-v1.0.1.def"
README = ROOT / "README.md"


class ContainerRecipeTests(unittest.TestCase):
    def test_recipe_installs_complete_pinned_kraken2_distribution(self):
        text = RECIPE.read_text(encoding="utf-8")

        self.assertRegex(text, r"KRAKEN2_VERSION=2\.1\.3\b")
        self.assertRegex(
            text,
            r"From: mambaorg/micromamba@sha256:[0-9a-f]{64}",
        )
        self.assertIn("/opt/kraken2", text)
        self.assertRegex(text, r"(?m)^\s*zlib1g-dev\b")
        self.assertRegex(text, r"(?m)^\s*centrifuge=1\.0\.4\.2\b")
        self.assertIn(
            "(cd /tmp/kraken2-src && ./install_kraken2.sh /opt/kraken2)",
            text,
        )
        self.assertRegex(text, r"export PATH=/opt/kraken2:[^\n]*\$PATH")
        self.assertNotRegex(text, r"(?m)^\s*cp\s+.*kraken2\*")

        for dependency in (
            "r-base",
            "r-inline",
            "r-gam",
            "r-rcpp",
            "r-rcppgsl",
            "r-ggplot2",
        ):
            self.assertRegex(text, rf"(?m)^\s*{re.escape(dependency)}(?:\s|\\|$)")
        self.assertIn("library(RcppGSL)", text)

        for executable in (
            "kraken2",
            "kraken2-build",
            "kraken2-inspect",
            "classify",
            "build_db",
            "k2",
            "estimate_capacity",
        ):
            self.assertRegex(
                text,
                rf"test -x /opt/kraken2/{re.escape(executable)}\b",
                msg=f"recipe does not verify {executable}",
            )

    def test_recipe_runs_functional_kraken2_smoke_test(self):
        text = RECIPE.read_text(encoding="utf-8")
        self.assertIn("/opt/sedimix/tests/smoke_test_kraken2.sh", text)

    def test_ci_builds_and_tests_the_versioned_recipe(self):
        workflow = (
            ROOT / ".github" / "workflows" / "container-smoke-test.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("container/sedimix-v1.0.1.def", workflow)
        self.assertIn("python3 tests/test_container_recipe.py -v", workflow)
        self.assertGreaterEqual(workflow.count('"tests/test_container_recipe.py"'), 2)
        self.assertGreaterEqual(workflow.count('"environment.yaml"'), 2)
        self.assertIn("apptainer build", workflow)
        self.assertIn("apptainer test", workflow)
        action_refs = re.findall(r"uses:\s+[^@\s]+@([^\s]+)", workflow)
        self.assertTrue(action_refs, "workflow must use at least one external action")
        for ref in action_refs:
            self.assertRegex(ref, r"^[0-9a-f]{40}$")

    def test_readme_links_peer_reviewed_article(self):
        text = README.read_text(encoding="utf-8")
        doi = "https://doi.org/10.1093/bioinformatics/btag004"
        self.assertIn(doi, text)
        self.assertNotIn("biorxiv.org/content/10.1101/2025.02.28.640818v1", text)
        self.assertIn("## Citation", text)
        self.assertIn("not yet published", text)
        self.assertIn(
            "library://jieruixu/sedimix/sedimix-v1.0:1.0.1",
            text,
        )
        self.assertNotIn("sedimix-v1.0.1:latest", text)


if __name__ == "__main__":
    unittest.main()
