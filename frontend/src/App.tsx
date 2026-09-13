import { useCallback, useEffect, useState } from "react";
import { predictVariant, PredictResponse } from "./api";

const AAS = "ARNDCQEGHILKMFPSTWYV".split("");

function getPath(): string {
  if (typeof window === "undefined") return "/";
  const p = window.location.pathname;
  if (p === "/predict" || p === "/predict/") return "/predict";
  if (p === "/methodology" || p === "/methodology/" || p === "/how-it-works") return "/methodology";
  return "/";
}

function navigate(to: string) {
  window.history.pushState({}, "", to);
  window.dispatchEvent(new PopStateEvent("popstate"));
  window.scrollTo(0, 0);
}

function useRoute(): string {
  const [path, setPath] = useState(getPath());
  useEffect(() => {
    const onPop = () => setPath(getPath());
    window.addEventListener("popstate", onPop);
    return () => window.removeEventListener("popstate", onPop);
  }, []);
  return path;
}

export default function App() {
  const route = useRoute();
  return (
    <div style={styles.page}>
      <style>{`
        * { box-sizing: border-box; }
        html, body { margin: 0; padding: 0; background: #f8fafc; }
        @media (max-width: 800px) {
          .vi-hero { grid-template-columns: 1fr !important; }
          .vi-features { grid-template-columns: 1fr !important; }
        }
        @media (max-width: 640px) {
          .vi-row { flex-direction: column !important; align-items: stretch !important; }
          .vi-row label { width: 100% !important; }
          .vi-navlinks { gap: 12px !important; font-size: 13px !important; }
        }
        @media (prefers-reduced-motion: reduce) {
          * { animation: none !important; transition: none !important; }
        }
        .skip-link:focus { left: 12px !important; top: 12px !important; width: auto !important; height: auto !important; overflow: visible !important; }
        a { color: #6d28d9; }
      `}</style>
      <a href="#main-content" className="skip-link" style={styles.skipLink}>Skip to main content</a>

      <nav style={styles.navbar}>
        <div style={styles.brand}>
          <button onClick={() => navigate("/")} style={styles.brandButton} aria-label="Variant Impact home">
            <span style={styles.brandMark}>VI</span>
            <span style={styles.brandName}>Variant Impact</span>
          </button>
        </div>
        <div className="vi-navlinks" style={styles.navLinks}>
          <NavLink label="Home" to="/" active={route === "/"} />
          <NavLink label="Check a variant" to="/predict" active={route === "/predict"} />
          <NavLink label="How it works" to="/methodology" active={route === "/methodology"} />
        </div>
      </nav>

      <main id="main-content">
        {route === "/predict" ? (
          <PredictPage />
        ) : route === "/methodology" ? (
          <MethodologyPage />
        ) : (
          <HomePage />
        )}
      </main>

      <footer style={styles.footer}>
        Variant Impact — for research and learning. Not medical advice. Always talk to a qualified
        clinician about health decisions.
      </footer>
    </div>
  );
}

function NavLink({ label, to, active }: { label: string; to: string; active: boolean }) {
  return (
    <button
      onClick={() => navigate(to)}
      aria-current={active ? "page" : undefined}
      style={active ? { ...styles.navLink, ...styles.navLinkActive } : styles.navLink}
    >
      {label}
    </button>
  );
}

/* ---------------- Home: hero + features + link only (no prediction form here) ---------------- */

function HomePage() {
  return (
    <div>
      <section className="vi-hero" style={styles.heroSection}>
        <div>
          <div style={styles.eyebrow}>Single-change gene scoring</div>
          <h1 style={styles.heroTitle}>
            Score a mutation, <em style={styles.heroEm}>see the reasoning.</em>
          </h1>
          <p style={styles.lede}>
            Variant Impact estimates how likely a single building-block swap in a protein is to
            affect how that protein works. It looks at the protein&apos;s shape and at how
            similar positions have changed across evolution — and every result comes with the
            reasons behind it in plain words.
          </p>
          <div style={styles.ctaRow}>
            <button onClick={() => navigate("/predict")} style={styles.buttonPrimary}>
              Check a variant
            </button>
            <button onClick={() => navigate("/methodology")} style={styles.buttonSecondary}>
              How it works
            </button>
          </div>
          <p style={styles.small}>Takes less than a minute. No account needed.</p>
        </div>
        <figure style={styles.heroVisual}>
          <img
            src="/hero.png"
            alt="Illustration of a folded protein beside a DNA double helix with one change highlighted"
            style={styles.heroImage}
          />
        </figure>
      </section>

      <section className="vi-features" style={styles.featureGrid} aria-label="Key features">
        <div style={styles.featureCard}>
          <span style={styles.featureIndex}>01</span>
          <h3 style={styles.featureTitle}>Looks at shape</h3>
          <p style={styles.featureText}>
            Reads the protein&apos;s real 3D shape to judge whether a position is buried inside
            or exposed on the surface.
          </p>
        </div>
        <div style={styles.featureCard}>
          <span style={styles.featureIndex}>02</span>
          <h3 style={styles.featureTitle}>Learns from evolution</h3>
          <p style={styles.featureText}>
            Weighs how unchanged a position has stayed across related proteins over long periods
            of time.
          </p>
        </div>
        <div style={styles.featureCard}>
          <span style={styles.featureIndex}>03</span>
          <h3 style={styles.featureTitle}>Explains every result</h3>
          <p style={styles.featureText}>
            Every score comes with the individual factors behind it, written so anyone can follow.
          </p>
        </div>
      </section>

      <section style={styles.homeCta} aria-label="Get started">
        <h2 style={styles.homeCtaTitle}>Have a change you want to understand?</h2>
        <p style={styles.homeCtaText}>
          Enter the protein, the position, and the swap — you will get a score plus a
          step-by-step explanation.
        </p>
        <button onClick={() => navigate("/predict")} style={styles.buttonPrimary}>
          Go to the checker
        </button>
      </section>
    </div>
  );
}

/* ---------------- Predict page ---------------- */

function PredictPage() {
  const [protein, setProtein] = useState("P53");
  const [position, setPosition] = useState("175");
  const [wt, setWt] = useState("R");
  const [mut, setMut] = useState("H");
  const [chain, setChain] = useState("A");
  const [result, setResult] = useState<PredictResponse | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  async function onPredict(e: React.FormEvent) {
    e.preventDefault();
    setError(null);
    setResult(null);
    setLoading(true);
    try {
      const pos = parseInt(position, 10);
      if (!Number.isFinite(pos) || pos < 1) throw new Error("Position must be a whole number of 1 or more.");
      const cleanProtein = protein.trim() || "UNKNOWN";
      if (!/^[A-Za-z0-9_\-]+$/.test(cleanProtein))
        throw new Error("Protein name may only use letters, numbers, dash or underscore.");
      const cleanChain = chain.trim().toUpperCase() || "A";
      if (!/^[A-Za-z0-9]$/.test(cleanChain)) throw new Error("Chain must be a single letter or number.");
      if (wt.toUpperCase() === mut.toUpperCase())
        throw new Error("The original and new building blocks must be different.");
      const res = await predictVariant({
        protein: cleanProtein,
        position: pos,
        wt: wt.toUpperCase(),
        mut: mut.toUpperCase(),
        chain: cleanChain,
      });
      setResult(res);
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : String(err);
      const status = (err as { status?: number })?.status;
      if (status === 503 || msg.toLowerCase().includes("not yet released")) {
        setError("Results are not available right now. Please try again in a little while.");
      } else {
        setError(friendlyError(msg));
      }
    } finally {
      setLoading(false);
    }
  }

  const plainLabel =
    result == null ? null : result.label === "damaging" ? "Likely to affect how the protein works" : "Unlikely to affect how the protein works";

  return (
    <div style={styles.narrowWrap}>
      <button onClick={() => navigate("/")} style={styles.backLink}>
        ← Back to home
      </button>
      <section style={styles.card} aria-labelledby="score-heading">
        <h1 id="score-heading" style={styles.cardTitleLarge}>Check a variant</h1>
        <p style={styles.cardDesc}>
          Enter a protein, a position, and the original and new building blocks. You will get a
          score plus the reasons behind it in plain words.
        </p>
        <form onSubmit={onPredict} style={styles.form} noValidate>
          <div style={styles.row} className="vi-row">
            <label htmlFor="protein-input" style={styles.label}>
              Gene or protein
              <input
                id="protein-input"
                name="protein"
                value={protein}
                onChange={(e) => setProtein(e.target.value)}
                style={styles.input}
                placeholder="P53"
                autoComplete="off"
              />
            </label>
            <label htmlFor="chain-input" style={styles.label}>
              Chain (usually A)
              <input
                id="chain-input"
                name="chain"
                value={chain}
                onChange={(e) => setChain(e.target.value)}
                style={{ ...styles.input, width: 90 }}
                maxLength={1}
                autoComplete="off"
              />
            </label>
            <label htmlFor="position-input" style={styles.label}>
              Position
              <input
                id="position-input"
                name="position"
                value={position}
                onChange={(e) => setPosition(e.target.value)}
                style={styles.input}
                type="number"
                min={1}
                step={1}
                inputMode="numeric"
              />
            </label>
          </div>
          <div style={styles.row} className="vi-row">
            <label htmlFor="wt-select" style={styles.label}>
              Original
              <select id="wt-select" name="wt" value={wt} onChange={(e) => setWt(e.target.value)} style={styles.input}>
                {AAS.map((a) => (
                  <option key={a} value={a}>{a}</option>
                ))}
              </select>
            </label>
            <span aria-hidden="true" style={styles.arrow}>→</span>
            <label htmlFor="mut-select" style={styles.label}>
              New
              <select id="mut-select" name="mut" value={mut} onChange={(e) => setMut(e.target.value)} style={styles.input}>
                {AAS.map((a) => (
                  <option key={a} value={a}>{a}</option>
                ))}
              </select>
            </label>
            <button type="submit" disabled={loading} style={styles.buttonPrimary}>
              {loading ? "Checking…" : "Check impact"}
            </button>
          </div>
          <div style={styles.example}>
            Example: protein P53, position 175, R → H. Most proteins use chain A.
          </div>
        </form>
        {error && <div style={styles.error} role="alert">{error}</div>}
      </section>

      {result && (
        <section style={styles.card} aria-live="polite" aria-label="Prediction result">
          <h2 style={styles.cardTitle}>Result for {result.variant}</h2>
          <div style={styles.scoreRow}>
            <div
              style={scoreCircleStyle(result.score)}
              role="img"
              aria-label={`Score ${result.score.toFixed(2)} out of 1, ${plainLabel}`}
            >
              {result.score.toFixed(2)}
            </div>
            <div style={{ flex: 1, minWidth: 220 }}>
              <div style={styles.resultLabel}>{plainLabel}</div>
              <div style={styles.resultHint}>
                Scores run from 0 to 1. Closer to 1 means more likely to affect how the protein works.
                Scores of 0.5 and above are grouped as likely to affect function.
              </div>
              <div
                style={styles.barWrap}
                role="progressbar"
                aria-valuenow={Math.round(result.score * 100)}
                aria-valuemin={0}
                aria-valuemax={100}
                aria-label="Likelihood score"
              >
                <div style={barFillStyle(result.score)} />
              </div>
            </div>
          </div>

          <h3 style={styles.subhead}>Why this score</h3>
          <PlainExplanations features={result.features} />
          <p style={styles.small}>
            This is a research estimate, not medical advice. Many real-world factors are not
            included. Want the full background?{" "}
            <button onClick={() => navigate("/methodology")} style={styles.inlineLink}>
              Read how it works
            </button>
            .
          </p>
        </section>
      )}
    </div>
  );
}

function friendlyError(msg: string): string {
  const lower = msg.toLowerCase();
  if (lower.includes("must be distinct") || lower.includes("synonymous"))
    return "The original and new building blocks must be different.";
  if (lower.includes("invalid amino acid"))
    return "Please choose valid building blocks from the lists.";
  if (lower.includes("protein must") || lower.includes("protein identifier"))
    return "Protein name may only use letters, numbers, dash or underscore.";
  if (lower.includes("chain must"))
    return "Chain must be a single letter or number (A is the usual choice).";
  if (lower.includes("position") || lower.includes("greater than or equal") || lower.includes("less than"))
    return "Position must be a whole number of 1 or more.";
  if (lower.includes("failed to fetch") || lower.includes("network") || lower.includes("health check failed"))
    return "Could not reach the scoring service. Check your connection and try again.";
  // Never leak file paths, revision hashes, or internal codes
  if (lower.includes("/") || lower.includes(".pkl") || lower.includes("revision") || lower.includes("checkpoint"))
    return "Something went wrong while scoring. Please try again.";
  return msg.length > 300 ? "Something went wrong while scoring. Please try again." : msg;
}

function PlainExplanations({ features }: { features: Record<string, number | string> }) {
  const num = (k: string): number | null => {
    const v = features[k];
    return typeof v === "number" && Number.isFinite(v) ? v : null;
  };
  const exposure = num("rsa");
  const conserved = num("conservation");
  const hydro = num("delta_hydrophobicity") ?? num("delta_hydro");
  const volume = num("delta_volume") ?? num("delta_vol");
  const grantham = num("grantham_distance") ?? num("grantham");
  const blosum = num("blosum62");
  const shape = features["secondary_structure"];

  const items: { title: string; text: string }[] = [];

  if (exposure != null)
    items.push({
      title: "How exposed the position is",
      text:
        exposure < 0.25
          ? `Buried deep inside the protein (${exposure.toFixed(2)} of 1). Changes in buried spots are more likely to disturb the structure.`
          : exposure < 0.6
          ? `Partly exposed (${exposure.toFixed(2)} of 1).`
          : `On or near the surface (${exposure.toFixed(2)} of 1). Surface changes are often easier for the protein to tolerate.`,
    });

  if (conserved != null)
    items.push({
      title: "How unchanged it stayed over evolution",
      text:
        conserved > 0.7
          ? `Very unchanged across related proteins (${conserved.toFixed(2)} of 1). When nature keeps a spot the same, a change there deserves attention.`
          : conserved > 0.4
          ? `Fairly unchanged across related proteins (${conserved.toFixed(2)} of 1).`
          : `Varies a lot across related proteins (${conserved.toFixed(2)} of 1). Spots that vary are often easier to change safely.`,
    });

  if (hydro != null)
    items.push({
      title: "Change in water preference",
      text:
        Math.abs(hydro) < 0.5
          ? `Almost the same (${hydro.toFixed(2)}). The swap keeps a similar liking for water.`
          : `Noticeably different (${hydro.toFixed(2)}). A big shift in water preference can disturb folding.`,
    });

  if (volume != null)
    items.push({
      title: "Change in size",
      text:
        Math.abs(volume) < 30
          ? `Similar in size (difference ${volume.toFixed(1)}). Fits in roughly the same space.`
          : `Quite different in size (difference ${volume.toFixed(1)}). A bigger size change is harder to fit.`,
    });

  if (grantham != null || blosum != null) {
    const bits: string[] = [];
    if (grantham != null)
      bits.push(grantham > 0.5 ? "chemically quite different" : grantham > 0.2 ? "somewhat different in chemistry" : "chemically similar");
    if (blosum != null)
      bits.push(blosum < 0 ? "a swap rarely seen in nature" : blosum === 0 ? "a swap seen from time to time" : "a swap often seen in nature");
    items.push({
      title: "How unusual the swap is",
      text: `The two building blocks are ${bits.join(" and ")}. Unusual swaps are more likely to matter.`,
    });
  }

  if (typeof shape === "string" && shape) {
    const shapePlain =
      shape === "H" ? "in a spiral section" : shape === "E" ? "in a flat sheet section" : "in a loop section";
    items.push({
      title: "Local shape",
      text: `The position sits ${shapePlain} of the protein. Shape helps decide how strict that spot is.`,
    });
  }

  if (items.length === 0) return <p style={styles.small}>No detailed breakdown is available for this result.</p>;

  return (
    <div style={styles.explainGrid}>
      {items.map((it) => (
        <div key={it.title} style={styles.explainCard}>
          <div style={styles.explainTitle}>{it.title}</div>
          <div style={styles.explainText}>{it.text}</div>
        </div>
      ))}
    </div>
  );
}

/* ---------------- Methodology page (technical detail lives here only) ---------------- */

function MethodologyPage() {
  const goPredict = useCallback(() => navigate("/predict"), []);
  return (
    <div style={styles.narrowWrap}>
      <button onClick={() => navigate("/")} style={styles.backLink}>← Back to home</button>
      <section style={styles.card} aria-labelledby="meth-heading">
        <h1 id="meth-heading" style={styles.cardTitleLarge}>How it works</h1>
        <p style={styles.cardDesc}>
          A short, honest tour of the data, the model, and the limits. This is the one place
          where we use technical names so experts can check our work.
        </p>

        <h2 style={styles.methHead}>What it does</h2>
        <p style={styles.methText}>
          For one missense change (one building-block swap at one position), the tool returns a
          score from 0 to 1. Closer to 1 means more likely to affect how the protein works.
          Scores of 0.5 and above are grouped as likely to affect function.
        </p>

        <h2 style={styles.methHead}>What it looks at</h2>
        <ul style={styles.list}>
          <li><strong>Protein shape:</strong> relative solvent accessibility and local shape (helix / sheet / loop) from real 3D structures, plus size and water-preference change between the two building blocks.</li>
          <li><strong>Evolutionary history:</strong> per-position conservation from alignments of related proteins (Shannon entropy / information content, 0 variable to 1 conserved).</li>
          <li><strong>Chemistry of the swap:</strong> Grantham distance and BLOSUM62 score, plus volume and hydrophobicity deltas (Kyte-Doolittle scale).</li>
          <li><strong>Combined score:</strong> a logistic-regression fusion model weighs all factors into one score.</li>
        </ul>

        <h2 style={styles.methHead}>Where the data comes from</h2>
        <ul style={styles.list}>
          <li><strong>Labels:</strong> ClinVar (NCBI) pathogenic / benign missense assertions (GRCh38 release).</li>
          <li><strong>Structures:</strong> AlphaFold DB predicted structures (e.g. AF-P04637 for human p53) and RCSB PDB entries; solvent access via DSSP (mkdssp) with a neighbor-count fallback.</li>
          <li><strong>Conservation:</strong> alignments of related (ortholog) sequences; per-column entropy turned into a 0–1 conservation score.</li>
        </ul>

        <h2 style={styles.methHead}>How well it works</h2>
        <p style={styles.methText}>
          In held-out testing the fusion model ranked a harmful change ahead of a harmless one
          about 8 times out of 10 (AUROC ≈ 0.80), ahead of conservation-only and
          chemistry-only baselines. That is a research-grade result, not a diagnosis.
        </p>

        <h2 style={styles.methHead}>Limits</h2>
        <ul style={styles.list}>
          <li>One change at a time; it does not judge combinations or whole-gene risk.</li>
          <li>It does not see clinical history, other genes, or environment.</li>
          <li>Shape and conservation estimates can be uncertain for poorly studied proteins.</li>
          <li>Never use this alone for health decisions — always talk to a qualified clinician.</li>
        </ul>

        <div style={styles.ctaRow}>
          <button onClick={goPredict} style={styles.buttonPrimary}>Try the checker</button>
        </div>
      </section>
    </div>
  );
}

/* ---------------- Light-theme styles ---------------- */

function scoreCircleStyle(s: number): React.CSSProperties {
  return {
    width: 76,
    height: 76,
    borderRadius: "50%",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    fontWeight: 800,
    fontSize: 17,
    color: "white",
    flexShrink: 0,
    background: s >= 0.5 ? (s > 0.75 ? "#991b1b" : "#dc2626") : s > 0.3 ? "#b45309" : "#15803d",
  };
}

function barFillStyle(s: number): React.CSSProperties {
  return {
    width: `${Math.round(s * 100)}%`,
    height: "100%",
    background: s >= 0.5 ? "#dc2626" : "#16a34a",
  };
}

const styles: Record<string, React.CSSProperties> = {
  page: {
    fontFamily: "'Inter', system-ui, -apple-system, sans-serif",
    background: "#f8fafc",
    minHeight: "100vh",
    color: "#0f172a",
  },
  skipLink: { position: "absolute", left: "-10000px", top: "auto", width: 1, height: 1, overflow: "hidden", background: "#0f172a", color: "white", padding: "8px 12px", borderRadius: 6, zIndex: 1000 } as React.CSSProperties,
  navbar: { maxWidth: 980, margin: "0 auto", padding: "16px 20px", display: "flex", justifyContent: "space-between", alignItems: "center", gap: 16, background: "#ffffff", borderBottom: "1px solid #e2e8f0", position: "sticky", top: 0, zIndex: 50 },
  brand: { display: "flex", alignItems: "center" },
  brandButton: { display: "flex", alignItems: "center", gap: 10, background: "none", border: "none", cursor: "pointer", padding: 0 },
  brandMark: { display: "inline-flex", alignItems: "center", justifyContent: "center", width: 34, height: 34, borderRadius: 9, background: "#6d28d9", color: "white", fontWeight: 700, fontSize: 13 },
  brandName: { fontWeight: 700, fontSize: 17, color: "#0f172a" },
  navLinks: { display: "flex", gap: 8, alignItems: "center" },
  navLink: { background: "none", border: "1px solid transparent", cursor: "pointer", fontSize: 14, fontWeight: 500, color: "#475569", padding: "8px 12px", borderRadius: 8 },
  navLinkActive: { color: "#6d28d9", background: "#f5f3ff", borderColor: "#ddd6fe", fontWeight: 700 },
  heroSection: { maxWidth: 980, margin: "0 auto", padding: "40px 20px", display: "grid", gridTemplateColumns: "1.1fr 0.9fr", gap: 36, alignItems: "center" },
  eyebrow: { textTransform: "uppercase", letterSpacing: "0.12em", fontSize: 12, fontWeight: 700, color: "#6d28d9", marginBottom: 12 },
  heroTitle: { fontWeight: 700, fontSize: 34, lineHeight: 1.15, margin: "0 0 16px", color: "#0f172a" },
  heroEm: { fontStyle: "italic", color: "#6d28d9" },
  lede: { color: "#475569", fontSize: 16, lineHeight: 1.6, maxWidth: 520, margin: "0 0 20px" },
  ctaRow: { display: "flex", gap: 12, flexWrap: "wrap", marginTop: 8 },
  heroVisual: { margin: 0 },
  heroImage: { width: "100%", maxWidth: 420, aspectRatio: "1 / 1", objectFit: "cover", borderRadius: 18, border: "1px solid #e2e8f0", boxShadow: "0 12px 32px rgba(15,23,42,0.10)", display: "block", background: "white" },
  featureGrid: { maxWidth: 980, margin: "0 auto", padding: "0 20px 8px", display: "grid", gridTemplateColumns: "repeat(3, 1fr)", gap: 16 },
  featureCard: { background: "#ffffff", border: "1px solid #e2e8f0", borderRadius: 14, padding: 22, boxShadow: "0 1px 2px rgba(15,23,42,0.05)" },
  featureIndex: { display: "block", fontSize: 13, fontWeight: 700, color: "#6d28d9", marginBottom: 8 },
  featureTitle: { fontWeight: 700, fontSize: 16, margin: "0 0 8px", color: "#0f172a" },
  featureText: { color: "#475569", fontSize: 14, lineHeight: 1.55, margin: 0 },
  homeCta: { maxWidth: 940, margin: "24px auto 32px", padding: "28px 24px", background: "#ffffff", border: "1px solid #e2e8f0", borderRadius: 14, textAlign: "center" },
  homeCtaTitle: { margin: "0 0 8px", fontSize: 20, color: "#0f172a" },
  homeCtaText: { margin: "0 0 16px", color: "#475569", fontSize: 14, lineHeight: 1.6 },
  narrowWrap: { maxWidth: 760, margin: "0 auto", padding: "28px 20px 32px" },
  backLink: { background: "none", border: "none", color: "#6d28d9", cursor: "pointer", fontSize: 14, fontWeight: 600, padding: "0 0 12px" },
  card: { background: "#ffffff", border: "1px solid #e2e8f0", borderRadius: 14, padding: 24, marginBottom: 18, boxShadow: "0 1px 2px rgba(15,23,42,0.05)" },
  cardTitle: { margin: "0 0 6px", fontSize: 18, color: "#0f172a", fontWeight: 700 },
  cardTitleLarge: { margin: "0 0 8px", fontSize: 24, color: "#0f172a", fontWeight: 800 },
  cardDesc: { margin: "0 0 16px", color: "#475569", fontSize: 14, lineHeight: 1.6 },
  form: { display: "flex", flexDirection: "column", gap: 14 },
  row: { display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-end" },
  label: { display: "flex", flexDirection: "column", gap: 6, fontSize: 13, fontWeight: 600, color: "#0f172a" },
  input: { padding: "10px 12px", borderRadius: 8, border: "1px solid #cbd5e1", fontSize: 14, minWidth: 100, color: "#0f172a", background: "#ffffff" },
  arrow: { alignSelf: "flex-end", paddingBottom: 10, fontSize: 20, color: "#64748b" },
  buttonPrimary: { background: "#6d28d9", color: "#ffffff", border: "none", padding: "11px 20px", borderRadius: 8, fontWeight: 700, fontSize: 14, cursor: "pointer" },
  buttonSecondary: { background: "#ffffff", color: "#6d28d9", border: "1px solid #ddd6fe", padding: "11px 20px", borderRadius: 8, fontWeight: 700, fontSize: 14, cursor: "pointer" },
  example: { fontSize: 12, color: "#64748b" },
  error: { background: "#fef2f2", border: "1px solid #fecaca", color: "#991b1b", padding: "10px 12px", borderRadius: 8, fontSize: 13, marginTop: 12 },
  scoreRow: { display: "flex", gap: 18, alignItems: "center", marginTop: 12, flexWrap: "wrap" },
  resultLabel: { fontSize: 18, fontWeight: 800, color: "#0f172a" },
  resultHint: { color: "#475569", fontSize: 13, lineHeight: 1.5, marginTop: 4 },
  subhead: { margin: "20px 0 10px", fontSize: 13, fontWeight: 700, textTransform: "uppercase", letterSpacing: 0.5, color: "#475569" },
  barWrap: { width: 220, maxWidth: "100%", height: 10, background: "#e2e8f0", borderRadius: 999, marginTop: 10, overflow: "hidden" },
  explainGrid: { display: "grid", gridTemplateColumns: "repeat(auto-fill, minmax(220px, 1fr))", gap: 12 },
  explainCard: { background: "#f8fafc", border: "1px solid #e2e8f0", borderRadius: 10, padding: "12px 14px" },
  explainTitle: { fontSize: 12, fontWeight: 700, color: "#6d28d9", textTransform: "uppercase", letterSpacing: 0.4 },
  explainText: { fontSize: 13, color: "#334155", lineHeight: 1.55, marginTop: 4 },
  methHead: { margin: "20px 0 8px", fontSize: 16, color: "#0f172a" },
  methText: { fontSize: 14, lineHeight: 1.65, color: "#334155", margin: "0 0 8px" },
  list: { fontSize: 14, lineHeight: 1.65, color: "#334155", paddingLeft: 20, margin: "0 0 8px" },
  small: { fontSize: 12, color: "#64748b", lineHeight: 1.5 },
  inlineLink: { background: "none", border: "none", padding: 0, color: "#6d28d9", fontWeight: 700, cursor: "pointer", fontSize: 12, textDecoration: "underline" },
  footer: { textAlign: "center", padding: "20px", fontSize: 12, color: "#64748b", borderTop: "1px solid #e2e8f0", background: "#ffffff", lineHeight: 1.5 },
};
