import { useEffect, useState } from "react";
import { fetchHealth, predictVariant, HealthResponse, PredictResponse } from "./api";

const AAS = "ARNDCQEGHILKMFPSTWYV".split("");

export default function App() {
  const [health, setHealth] = useState<HealthResponse | null>(null);
  const [healthError, setHealthError] = useState<string | null>(null);
  const [protein, setProtein] = useState("P53");
  const [position, setPosition] = useState("175");
  const [wt, setWt] = useState("R");
  const [mut, setMut] = useState("H");
  const [chain, setChain] = useState("A");
  const [result, setResult] = useState<PredictResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  useEffect(() => {
    fetchHealth()
      .then(setHealth)
      .catch((e) => setHealthError(e.message));
  }, []);

  async function onPredict(e: React.FormEvent) {
    e.preventDefault();
    setError(null);
    setResult(null);
    setLoading(true);
    try {
      const pos = parseInt(position, 10);
      if (!Number.isFinite(pos) || pos < 1) throw new Error("Position must be a positive integer (≥1)");
      const cleanProtein = protein.trim() || "UNKNOWN";
      if (!/^[A-Za-z0-9_\-]+$/.test(cleanProtein)) throw new Error("Protein identifier must be alphanumeric, underscore or dash");
      const cleanChain = chain.trim().toUpperCase() || "A";
      if (!/^[A-Za-z0-9]$/.test(cleanChain)) throw new Error("Chain must be a single alphanumeric character");
      const res = await predictVariant({ protein: cleanProtein, position: pos, wt: wt.toUpperCase(), mut: mut.toUpperCase(), chain: cleanChain });
      setResult(res);
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : String(err);
      // Enhance 503 abstention messaging
      const status = (err as { status?: number })?.status;
      if (status === 503 || msg.toLowerCase().includes("not yet released")) {
        setError("Model not yet released — abstaining. The backend release gate is closed; predictions are disabled until an approved artifact is released.");
      } else {
        setError(msg);
      }
    } finally {
      setLoading(false);
    }
  }

  const modelNotReleased = health ? !health.model_loaded : false;

  return (
    <div style={styles.page}>
      <style>{`
        @media (max-width: 640px) {
          .vi-row { flex-direction: column !important; align-items: stretch !important; }
          .vi-row label { width: 100% !important; }
          .vi-health { width: 100%; }
        }
        @media (prefers-reduced-motion: reduce) {
          * { animation: none !important; transition: none !important; }
        }
        .skip-link:focus { left: 12px !important; top: 12px !important; width: auto !important; height: auto !important; overflow: visible !important; }
      `}</style>
      <a href="#main-content" className="skip-link" style={styles.skipLink}>Skip to main content</a>
      <header style={styles.header} role="banner">
        <div style={styles.headerInner}>
          <div>
            <h1 style={styles.title}>Variant Impact</h1>
            <p style={styles.subtitle}>Structure/evolution-aware missense scoring — per-variant damaging prediction</p>
          </div>
          <div style={styles.healthBadge} role="status" aria-live="polite" aria-label="Model health status" className="vi-health">
            {healthError ? (
              <span style={{ color: "#991b1b" }}>Health: unreachable — backend not reachable</span>
            ) : health ? (
              <span style={{ color: health.model_loaded ? "#166534" : "#92400e" }}>
                Model: {health.model_loaded ? `loaded (rev ${health.model_revision})` : "not yet released"}
              </span>
            ) : (
              <span>Checking model…</span>
            )}
          </div>
        </div>
      </header>

      {modelNotReleased && (
        <div style={styles.banner} role="alert" aria-live="polite">
          <strong>Model not yet released — abstaining.</strong> The backend release gate is closed
          (<code>MODEL_RELEASE_APPROVED</code> / <code>APPROVED_ARTIFACT_REVISION</code> not set). Predictions will return 503 until an approved artifact is released. This is honest abstention, not a silent fallback.
        </div>
      )}

      <main id="main-content" style={styles.main}>
        <section style={styles.card} aria-labelledby="score-heading">
          <h2 id="score-heading" style={styles.cardTitle}>Score a missense variant</h2>
          <p style={styles.cardDesc}>
            Enter protein, position (1-indexed), wild-type and mutant residues. The model combines solvent accessibility (RSA),
            conservation (MSA entropy), and physicochemical deltas into a damaging score.
          </p>
          <form onSubmit={onPredict} style={styles.form} noValidate aria-describedby="form-help">
            <p id="form-help" style={{ position: "absolute", left: "-10000px", width: 1, height: 1, overflow: "hidden" }}>
              All fields are required. Position must be positive integer. WT and Mutant must be distinct amino acids.
            </p>
            <div style={styles.row} className="vi-row">
              <label htmlFor="protein-input" style={styles.label}>
                Protein / Gene
                <input id="protein-input" name="protein" value={protein} onChange={(e) => setProtein(e.target.value)} style={styles.input} placeholder="P53" aria-required="true" autoComplete="off" />
              </label>
              <label htmlFor="chain-input" style={styles.label}>
                Chain
                <input id="chain-input" name="chain" value={chain} onChange={(e) => setChain(e.target.value)} style={{ ...styles.input, width: 60 }} maxLength={1} aria-label="Chain ID single letter" aria-required="true" />
              </label>
              <label htmlFor="position-input" style={styles.label}>
                Position
                <input id="position-input" name="position" value={position} onChange={(e) => setPosition(e.target.value)} style={styles.input} type="number" min={1} step={1} aria-required="true" inputMode="numeric" />
              </label>
            </div>
            <div style={styles.row} className="vi-row">
              <label htmlFor="wt-select" style={styles.label}>
                WT
                <select id="wt-select" name="wt" value={wt} onChange={(e) => setWt(e.target.value)} style={styles.input} aria-label="Wild-type amino acid">
                  {AAS.map((a) => (
                    <option key={a} value={a}>
                      {a}
                    </option>
                  ))}
                </select>
              </label>
              <span aria-hidden="true" style={{ alignSelf: "flex-end", paddingBottom: 10, fontSize: 20 }}>→</span>
              <label htmlFor="mut-select" style={styles.label}>
                Mutant
                <select id="mut-select" name="mut" value={mut} onChange={(e) => setMut(e.target.value)} style={styles.input} aria-label="Mutant amino acid">
                  {AAS.map((a) => (
                    <option key={a} value={a}>
                      {a}
                    </option>
                  ))}
                </select>
              </label>
              <button type="submit" disabled={loading} style={styles.button} aria-busy={loading} aria-label={loading ? "Scoring variant, please wait" : "Predict impact of variant"}>
                {loading ? "Scoring…" : "Predict impact"}
              </button>
            </div>
            <div style={styles.example}>
              Example: <code>P53_R175H</code> — protein P53, position 175, R→H
            </div>
          </form>
          {error && <div style={styles.error} role="alert" aria-live="assertive">{error}</div>}
          <div aria-live="polite" aria-atomic="true" style={{ position: "absolute", left: "-10000px" }}>{loading ? "Scoring in progress" : ""}</div>
        </section>

        {result && (
          <section style={styles.card} aria-live="polite" aria-atomic="true" aria-label="Prediction result">
            <h3 style={styles.cardTitle}>Result: {result.variant}</h3>
            <div style={styles.scoreRow}>
              <div style={scoreCircleStyle(result.score)} role="img" aria-label={`Score ${result.score.toFixed(3)} out of 1, ${result.label}`}>{result.score.toFixed(3)}</div>
              <div>
                <div style={{ fontSize: 20, fontWeight: 700, color: result.label === "damaging" ? "#991b1b" : "#166534" }}>
                  {result.label.toUpperCase()}
                </div>
                <div style={{ color: "#475569", fontSize: 13 }}>
                  Score in [0,1]; threshold 0.5. Model rev: {result.model_revision}
                </div>
                <div style={styles.barWrap} role="progressbar" aria-valuenow={Math.round(result.score*100)} aria-valuemin={0} aria-valuemax={100} aria-label="Damaging score">
                  <div style={barFillStyle(result.score)} />
                </div>
              </div>
            </div>
            <h4 style={styles.subhead}>Per-feature explanation</h4>
            <div style={styles.grid}>
              <Feature label="RSA" value={result.features.rsa} hint="0 buried → 1 exposed" />
              <Feature label="Conservation" value={result.features.conservation} hint="0 variable → 1 conserved" />
              <Feature label="Δhydrophobicity" value={result.features.delta_hydrophobicity} hint="Kyte-Doolittle" />
              <Feature label="Δvolume" value={result.features.delta_volume} hint="Å³" />
              <Feature label="Grantham distance" value={result.features.grantham_distance} hint="0 similar → >100 dissimilar" />
              <Feature label="BLOSUM62" value={result.features.blosum62} hint="positive = common substitution" />
              <Feature label="Secondary structure" value={result.features.secondary_structure as string} hint="H/E/C" />
            </div>
            {result.features.note && <p style={styles.note}>{String(result.features.note)}</p>}
          </section>
        )}

        <section style={styles.card}>
          <h3 style={styles.cardTitle}>How it works</h3>
          <ul style={styles.list}>
            <li>
              <strong>Structural context:</strong> RSA and secondary structure from DSSP (or geometric neighbor-count proxy in sandbox) via{" "}
              <code>Bio.PDB</code> on real PDB/AlphaFold structures.
            </li>
            <li>
              <strong>Evolutionary conservation:</strong> Shannon-entropy information content from a real MSA FASTA.
            </li>
            <li>
              <strong>Physicochemical deltas:</strong> real Kyte-Doolittle, vdW volumes, Grantham distance, BLOSUM62.
            </li>
            <li>
              <strong>Fusion:</strong> logistic regression combining all features; evaluated AUROC/AUPRC vs single-feature baselines.
            </li>
            <li>
              <strong>Real-run:</strong> <code>python -m data_pipeline.cli --structure-path &lt;file&gt; --msa-path &lt;file&gt; --predict --variant P53_R175H</code>
            </li>
          </ul>
          <p style={styles.small}>
            Data sources: RCSB PDB, AlphaFold DB, DSSP, ClinVar, phyloP/phastCons — see <code>docs/data_sources.md</code>.
          </p>
        </section>
      </main>

      <footer style={styles.footer}>Variant Impact — research-grade, not clinical advice.</footer>
    </div>
  );
}

function Feature({ label, value, hint }: { label: string; value: unknown; hint?: string }) {
  const display = typeof value === "number" ? value.toFixed(3) : String(value ?? "—");
  return (
    <div style={styles.featureCell}>
      <div style={styles.featureLabel}>{label}</div>
      <div style={styles.featureValue}>{display}</div>
      {hint && <div style={styles.featureHint}>{hint}</div>}
    </div>
  );
}

const styles: Record<string, React.CSSProperties> = {
  page: { fontFamily: "Inter, system-ui, -apple-system, sans-serif", background: "#f8fafc", minHeight: "100vh", color: "#0f172a" },
  header: { background: "white", borderBottom: "1px solid #e2e8f0", padding: "18px 24px" },
  headerInner: { maxWidth: 980, margin: "0 auto", display: "flex", justifyContent: "space-between", alignItems: "center", gap: 16, flexWrap: "wrap" },
  title: { margin: 0, fontSize: 22, letterSpacing: -0.5 },
  subtitle: { margin: "4px 0 0", color: "#475569", fontSize: 13 },
  healthBadge: { fontSize: 13, background: "#f1f5f9", padding: "6px 12px", borderRadius: 999, border: "1px solid #e2e8f0" },
  banner: { maxWidth: 980, margin: "16px auto 0", background: "#fef3c7", border: "1px solid #fcd34d", padding: "12px 16px", borderRadius: 8, fontSize: 13, color: "#78350f" },
  main: { maxWidth: 980, margin: "20px auto", padding: "0 16px", display: "flex", flexDirection: "column", gap: 18 },
  card: { background: "white", border: "1px solid #e2e8f0", borderRadius: 12, padding: 20, boxShadow: "0 1px 2px rgba(0,0,0,0.04)" },
  cardTitle: { margin: "0 0 6px", fontSize: 16, color: "#0f172a" },
  cardDesc: { margin: "0 0 14px", color: "#475569", fontSize: 13, lineHeight: 1.5 },
  form: { display: "flex", flexDirection: "column", gap: 12 },
  row: { display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-end" },
  label: { display: "flex", flexDirection: "column", gap: 4, fontSize: 12, fontWeight: 600, color: "#1e293b" },
  input: { padding: "8px 10px", borderRadius: 8, border: "1px solid #94a3b8", fontSize: 14, minWidth: 80, color: "#0f172a", background: "white" },
  button: { background: "#0f172a", color: "white", border: "none", padding: "10px 18px", borderRadius: 8, fontWeight: 600, cursor: "pointer" },
  example: { fontSize: 12, color: "#64748b" },
  error: { background: "#fef2f2", border: "1px solid #fecaca", color: "#7f1d1d", padding: "10px 12px", borderRadius: 8, fontSize: 13 },
  scoreRow: { display: "flex", gap: 18, alignItems: "center", marginTop: 8, flexWrap: "wrap" },
  subhead: { margin: "18px 0 8px", fontSize: 13, textTransform: "uppercase", letterSpacing: 0.5, color: "#334155" },
  grid: { display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(150px, 1fr))", gap: 10 },
  featureCell: { background: "#f8fafc", border: "1px solid #e2e8f0", borderRadius: 8, padding: "10px 12px" },
  featureLabel: { fontSize: 11, fontWeight: 700, color: "#475569", textTransform: "uppercase", letterSpacing: 0.4 },
  featureValue: { fontSize: 18, fontWeight: 700, marginTop: 2, color: "#0f172a" },
  featureHint: { fontSize: 11, color: "#64748b", marginTop: 2 },
  barWrap: { width: 200, height: 8, background: "#e2e8f0", borderRadius: 999, marginTop: 6, overflow: "hidden" },
  note: { fontSize: 12, color: "#475569", fontStyle: "italic", marginTop: 10 },
  list: { fontSize: 13, lineHeight: 1.6, color: "#334155", paddingLeft: 18 },
  small: { fontSize: 12, color: "#64748b" },
  footer: { textAlign: "center", padding: "18px 0 28px", fontSize: 12, color: "#64748b" },
  skipLink: { position: "absolute", left: "-10000px", top: "auto", width: 1, height: 1, overflow: "hidden", background: "#0f172a", color: "white", padding: "8px 12px", borderRadius: 6, zIndex: 1000 } as React.CSSProperties,
};

function scoreCircleStyle(s: number): React.CSSProperties {
  return {
    width: 72,
    height: 72,
    borderRadius: "50%",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    fontWeight: 800,
    fontSize: 16,
    color: "white",
    background: s >= 0.5 ? (s > 0.75 ? "#991b1b" : "#dc2626") : s > 0.3 ? "#a16207" : "#15803d",
  };
}

function barFillStyle(s: number): React.CSSProperties {
  return {
    width: `${Math.round(s * 100)}%`,
    height: "100%",
    background: s >= 0.5 ? "#dc2626" : "#16a34a",
  };
}
