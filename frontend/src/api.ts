export interface HealthResponse {
  status: string;
  model_loaded: boolean;
  model_revision?: string | null;
}

export interface PredictResponse {
  variant: string;
  score: number;
  label: string;
  features: Record<string, number | string>;
  model_revision: string;
}

// Default backend host (Hugging Face Space). Env vars override at build time.
const DEFAULT_API = "https://bhumika-tewari-282006-variant-functional-impact-ai-api.hf.space";

// Vite exposes env via import.meta.env; VITE_API_URL is preferred, VITE_API_BASE is deprecated alias
const _env = (import.meta as unknown as { env: Record<string, string | undefined> }).env;
const BASE = (_env.VITE_API_URL || _env.VITE_API_BASE || DEFAULT_API).replace(/\/$/, "");

export async function fetchHealth(): Promise<HealthResponse> {
  const res = await fetch(`${BASE}/health`);
  if (!res.ok) throw new Error(`Health check failed: ${res.status}`);
  return res.json();
}

export async function predictVariant(payload: {
  protein: string;
  position: number;
  wt: string;
  mut: string;
  chain?: string;
}): Promise<PredictResponse> {
  const res = await fetch(`${BASE}/predict`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  if (!res.ok) {
    const body = await res.json().catch(() => ({ detail: res.statusText }));
    let msg: string;
    if (Array.isArray(body.detail)) {
      // Pydantic validation error list
      msg = body.detail.map((d: { msg?: string; loc?: unknown }) => d.msg || JSON.stringify(d)).join("; ");
    } else if (typeof body.detail === "string") {
      msg = body.detail;
    } else if (body.detail) {
      msg = JSON.stringify(body.detail);
    } else {
      msg = `Predict failed: ${res.status} ${res.statusText}`;
    }
    // Include status for caller to distinguish abstention
    const err = new Error(msg) as Error & { status?: number };
    err.status = res.status;
    throw err;
  }
  return res.json();
}
