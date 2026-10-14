scpPower 1.0.0 provides calibration, hierarchical simulation and power estimates for two-condition single-cell proteomics, with nominal instrument hours. The release preserves the registered scientific decisions and accepted M2, M4 and M4b verdicts.

M3b regenerates the 630 Specht and Leduc points with the amendment-I endpoint masks, original designs and seeds, and 200 replicates per point. The 315 non-TMT grid points are reused without rerunning. No cheapest-design row changed, and no design entered or left the set reaching 0.8. The full comparison and all refreshed tables, plots and Reading the grid tables are included in the unified evidence attachment.

The README contains worked examples for carrier TMT, plexDIA and label-free designs. The vignette is built and checked in CI; its displayed numerical results are computed at build time. The pilot protocol describes donor variance, ratio standards and blank controls. No registered quantity or validation gate changed.

Real-cell carrier-TMT validation reaches observed power of about 0.35. Donor-variance transfer and compression remain assumptions; cell-depth dependence is omitted. Q Exactive and timsTOF lack sufficient public cells for validation at useful power. Read the full limitations before using design estimates.

### Zenodo: three author steps

1. Sign in at [Zenodo](https://zenodo.org) with your GitHub account and connect that account if needed.
2. Open the profile menu → GitHub, click **Sync now**, and switch on **yhengaa/scpPower**. [Official repository-enabling instructions](https://help.zenodo.org/docs/github/enable-repository/).
3. Create the GitHub release after enabling. If v1.0.0 was already released, re-tag v1.0.0 at the **same green checked commit** and recreate the release after enabling, preserving its evidence attachments; wait for Zenodo to process it, then copy its DOI. Add the DOI badge and the DOI to CITATION.cff in a v1.0.1 commit. [Official release-archiving instructions](https://help.zenodo.org/docs/github/archive-software/github-upload/).

Zenodo integration is an author action and was not enabled by the assistant. No DOI exists in this bundle; no DOI or badge is invented. The tag is created only after R CMD check succeeds at its exact commit.

Code was written with the assistance of GPT (OpenAI) and independently reviewed with Claude (Anthropic), both under the author's direction. All scientific decisions, registered criteria, data choices and verifications are the author's.
