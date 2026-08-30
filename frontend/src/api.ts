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

const BASE = import.meta.env.VITE_API_BASE || "";

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
    throw new Error(body.detail || `Predict failed: ${res.status}`);
  }
  return res.json();
}
