# Production Deployment — Fast Startup with Sysimage

Julia's JIT compilation causes slow startup (~30-60s for EvolvingDomains + Gridap). 
A **sysimage** pre-compiles everything once, giving near-instant startup afterward.

---

## Quick Start

### 1. Build the Sysimage (One Time)

```bash
cd /path/to/EvolvingDomains
julia --project=. scripts/build_sysimage.jl
```

**Time:** 5-20 minutes (depends on machine)  
**Output:** `evolving_domains.so` (Linux/Mac) or `evolving_domains.dll` (Windows)

### 2. Run Simulations (Every Time)

```bash
julia --sysimage=evolving_domains.so your_script.jl
```

Startup will now take **< 2 seconds** instead of 30-60s.

---

## CI/CD Integration

### Docker

```dockerfile
FROM julia:1.10

WORKDIR /app
COPY . .

# Build sysimage during image build
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()'
RUN julia --project=. scripts/build_sysimage.jl

# Use sysimage for all runs
ENTRYPOINT ["julia", "--sysimage=evolving_domains.so"]
```

### GitHub Actions

```yaml
- name: Build Sysimage
  run: julia --project=. scripts/build_sysimage.jl

- name: Cache Sysimage
  uses: actions/cache@v3
  with:
    path: evolving_domains.so
    key: sysimage-${{ hashFiles('Project.toml', 'Manifest.toml') }}
```

---

## When to Rebuild

Rebuild the sysimage when:
- `Project.toml` or `Manifest.toml` changes (package updates)
- Julia version changes
- You add new workflows to `warmup_sysimage.jl`

---

## Files

| File | Purpose |
|------|---------|
| `build_sysimage.jl` | Main build script |
| `warmup_sysimage.jl` | Exercises code paths to compile (edit to add your workflows) |

---

## Troubleshooting

**"Method not found" or still slow for some operations?**  
→ Add those operations to `warmup_sysimage.jl` and rebuild.

**Sysimage too large?**  
→ Remove unused packages from the `PACKAGES` list in `build_sysimage.jl`.

**Crashes on different machine?**  
→ Change `cpu_target = "generic"` to match target architecture, or build on deployment machine.
