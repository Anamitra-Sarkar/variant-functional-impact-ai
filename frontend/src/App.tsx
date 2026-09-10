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
        setError("This prediction isn't available yet — our team is finishing validation before enabling live results.");
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
        @media (max-width: 800px) {
          .vi-hero { grid-template-columns: 1fr !important; }
          .vi-features { grid-template-columns: 1fr !important; }
        }
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

      <nav style={styles.navbar}>
        <div style={styles.brand}>
          <span style={styles.brandMark}>VI</span>
          <span style={styles.brandName}>Variant Impact</span>
        </div>
        <div role="status" aria-live="polite" aria-label="Model health status" className="vi-health" style={styles.healthBadge}>
          {healthError ? (
            <span style={{ color: "#f87171" }}>Not available</span>
          ) : health ? (
            <span style={{ color: health.model_loaded ? "#6fcf97" : "#c4bce0" }}>
              {health.model_loaded ? "Live predictions" : "Preview mode"}
            </span>
          ) : (
            <span>Checking…</span>
          )}
        </div>
      </nav>

      <section className="vi-hero" style={styles.heroSection}>
        <div style={styles.heroCopy}>
          <div style={styles.eyebrow}>Missense variant scoring</div>
          <h1 style={styles.heroTitle}>
            Score a mutation, <em style={styles.heroEm}>see the reasoning.</em>
          </h1>
          <p style={styles.lede}>
            Variant Impact estimates how damaging a single amino-acid change is likely to be, by combining a
            protein's own structure with how conserved that position is across evolution — every score comes with
            the evidence behind it.
          </p>
        </div>
        <figure style={styles.heroVisual}>
          <img
            src="/hero.png"
            alt="3D ribbon illustration of a folded protein with beta sheets and alpha helices in purple beside a DNA double helix with two base pairs highlighted in glowing pink to indicate a missense variant"
            style={{ width: "100%", aspectRatio: "1 / 1", objectFit: "cover", borderRadius: 18, border: "1px solid rgba(255,255,255,0.12)", boxShadow: "0 24px 60px rgba(0,0,0,0.5)", display: "block" }}
          />
        </figure>
      </section>

      <section className="vi-features" style={styles.featureGrid}>
        <div style={styles.featureCard}>
          <span style={styles.featureIndex}>01</span>
          <h3 style={styles.featureTitle}>Structure-aware</h3>
          <p style={styles.featureText}>Reads real protein structure to judge whether a position is buried or exposed.</p>
        </div>
        <div style={styles.featureCard}>
          <span style={styles.featureIndex}>02</span>
          <h3 style={styles.featureTitle}>Evolution-informed</h3>
          <p style={styles.featureText}>Weighs how conserved a position is across related proteins throughout evolution.</p>
        </div>
        <div style={styles.featureCard}>
          <span style={styles.featureIndex}>03</span>
          <h3 style={styles.featureTitle}>Fully explained</h3>
          <p style={styles.featureText}>Every score comes with the individual factors that contributed to it.</p>
        </div>
      </section>

      {modelNotReleased && (
        <div style={styles.banner} role="alert" aria-live="polite">
          <strong>Predictions aren't available yet.</strong> Our team is finishing validation before enabling live results.
        </div>
      )}

      <main id="main-content" style={styles.main}>
        <section style={styles.card} aria-labelledby="score-heading">
          <h2 id="score-heading" style={styles.cardTitle}>Score a missense variant</h2>
          <p style={styles.cardDesc}>
            Enter a protein, position, and the original and new amino acid to see how damaging that change is likely to be, and why.
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
              <strong>Structural context:</strong> how buried or exposed the position is within the protein's real 3D structure.
            </li>
            <li>
              <strong>Evolutionary conservation:</strong> how consistently that position is preserved across related proteins.
            </li>
            <li>
              <strong>Physicochemical change:</strong> how different the original and new amino acids are in size, charge, and chemistry.
            </li>
            <li>
              <strong>Combined scoring:</strong> all factors are weighed together into a single damaging-likelihood score.
            </li>
          </ul>
          <p style={styles.small}>
            Built on real protein structure and evolutionary conservation data from public research databases.
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
  page: {
    fontFamily: "'Inter', system-ui, -apple-system, sans-serif",
    background:
      "linear-gradient(rgba(14,11,33,0.92), rgba(14,11,33,0.97)), url('/hero.png')",
    backgroundSize: "cover",
    backgroundPosition: "top center",
    backgroundAttachment: "fixed",
    minHeight: "100vh",
    color: "#f1eefc",
  },
  navbar: { maxWidth: 980, margin: "0 auto", padding: "24px 20px 0", display: "flex", justifyContent: "space-between", alignItems: "center" },
  brand: { display: "flex", alignItems: "center", gap: 10 },
  brandMark: { display: "inline-flex", alignItems: "center", justifyContent: "center", width: 34, height: 34, borderRadius: 9, background: "linear-gradient(135deg, #a78bfa, #4c1d95)", color: "white", fontWeight: 700, fontSize: 13 },
  brandName: { fontFamily: "'Fraunces', serif", fontWeight: 600, fontSize: 18 },
  heroSection: { maxWidth: 980, margin: "0 auto", padding: "32px 20px 40px", display: "grid", gridTemplateColumns: "1.1fr 0.9fr", gap: 40, alignItems: "center" },
  heroCopy: {},
  eyebrow: { textTransform: "uppercase", letterSpacing: "0.14em", fontSize: 12, fontWeight: 700, color: "#a78bfa", marginBottom: 14 },
  heroTitle: { fontFamily: "'Fraunces', serif", fontWeight: 600, fontSize: 32, lineHeight: 1.15, margin: "0 0 18px" },
  heroEm: { fontStyle: "italic", color: "#a78bfa" },
  lede: { color: "#c4bce0", fontSize: 16, lineHeight: 1.6, maxWidth: 480, margin: 0 },
  heroVisual: { margin: 0 },
  featureGrid: { maxWidth: 980, margin: "0 auto 32px", padding: "0 20px", display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 18 },
  featureCard: { background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.12)", borderRadius: 14, padding: 22 },
  featureIndex: { display: "block", fontFamily: "'Fraunces', serif", fontSize: 13, color: "#a78bfa", marginBottom: 8 },
  featureTitle: { fontFamily: "'Fraunces', serif", fontWeight: 600, fontSize: 16, margin: "0 0 8px" },
  featureText: { color: "#c4bce0", fontSize: 13, lineHeight: 1.55, margin: 0 },
  title: { margin: 0, fontSize: 22, letterSpacing: -0.5 },
  subtitle: { margin: "4px 0 0", color: "#475569", fontSize: 13 },
  healthBadge: { fontSize: 13, background: "rgba(255,255,255,0.06)", padding: "6px 12px", borderRadius: 999, border: "1px solid rgba(255,255,255,0.12)" },
  banner: { maxWidth: 980, margin: "0 auto 20px", background: "rgba(217,175,86,0.1)", border: "1px solid rgba(217,175,86,0.5)", padding: "12px 16px", borderRadius: 8, fontSize: 13, color: "#f1eefc" },
  hero: { maxWidth: 980, margin: "0 auto", padding: "18px 16px 0" },
  heroImage: { width: "100%", maxHeight: 320, objectFit: "contain", display: "block", borderRadius: 12, background: "white", border: "1px solid #e2e8f0" },
  main: { maxWidth: 980, margin: "0 auto 20px", padding: "0 20px", display: "flex", flexDirection: "column", gap: 18 },
  card: { background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.12)", borderRadius: 12, padding: 20 },
  cardTitle: { margin: "0 0 6px", fontSize: 16, color: "#f1eefc", fontFamily: "'Fraunces', serif", fontWeight: 600 },
  cardDesc: { margin: "0 0 14px", color: "#c4bce0", fontSize: 13, lineHeight: 1.5 },
  form: { display: "flex", flexDirection: "column", gap: 12 },
  row: { display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-end" },
  label: { display: "flex", flexDirection: "column", gap: 4, fontSize: 12, fontWeight: 600, color: "#f1eefc" },
  input: { padding: "8px 10px", borderRadius: 8, border: "1px solid rgba(255,255,255,0.2)", fontSize: 14, minWidth: 80, color: "#f1eefc", background: "rgba(0,0,0,0.2)" },
  button: { background: "#a78bfa", color: "#1a1332", border: "none", padding: "10px 18px", borderRadius: 8, fontWeight: 600, cursor: "pointer" },
  example: { fontSize: 12, color: "#8b81ab" },
  error: { background: "rgba(248,113,113,0.1)", border: "1px solid rgba(248,113,113,0.4)", color: "#fecaca", padding: "10px 12px", borderRadius: 8, fontSize: 13 },
  scoreRow: { display: "flex", gap: 18, alignItems: "center", marginTop: 8, flexWrap: "wrap" },
  subhead: { margin: "18px 0 8px", fontSize: 13, textTransform: "uppercase", letterSpacing: 0.5, color: "#334155" },
  grid: { display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(150px, 1fr))", gap: 10 },
  featureCell: { background: "rgba(0,0,0,0.2)", border: "1px solid rgba(255,255,255,0.12)", borderRadius: 8, padding: "10px 12px" },
  featureLabel: { fontSize: 11, fontWeight: 700, color: "#c4bce0", textTransform: "uppercase", letterSpacing: 0.4 },
  featureValue: { fontSize: 18, fontWeight: 700, marginTop: 2, color: "#f1eefc" },
  featureHint: { fontSize: 11, color: "#8b81ab", marginTop: 2 },
  barWrap: { width: 200, height: 8, background: "#e2e8f0", borderRadius: 999, marginTop: 6, overflow: "hidden" },
  note: { fontSize: 12, color: "#c4bce0", fontStyle: "italic", marginTop: 10 },
  list: { fontSize: 13, lineHeight: 1.6, color: "#c4bce0", paddingLeft: 18 },
  small: { fontSize: 12, color: "#8b81ab" },
  footer: { textAlign: "center", padding: "18px 0 28px", fontSize: 12, color: "#8b81ab" },
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
