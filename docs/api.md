# HTTP service

`bellhop-server` is a stateless, synchronous HTTP facade. Each simulation
request supplies one [self-contained JSON case](json-input.md), waits for the
solver, and receives the completed HDF5 file directly. The service does not
retain cases or results.

## Start the service

```console
RUST_LOG=info cargo run -p bellhop-server --release
```

The default listener is `0.0.0.0:8080`. Generated documentation is public:

- OpenAPI 3 document: `GET /openapi.json`
- Swagger UI: `GET /docs`
- Readiness: `GET /healthz`

The service provides plain HTTP only. Terminate TLS in a reverse proxy or load
balancer when traffic crosses a trusted-network boundary.

## Endpoints

### `POST /v1/validate`

Accepts `application/json`. A valid case returns:

```json
{
  "valid": true,
  "warnings": []
}
```

Validation parses the complete document and checks its semantic consistency;
it does not run the numerical solver.

### `POST /v1/run`

Accepts the same JSON document. Success returns the HDF5 bytes with:

```text
Content-Type: application/x-hdf5
Content-Disposition: attachment; filename="bellhop-result.h5"
```

Example:

```console
curl --fail-with-body \
  -H 'Content-Type: application/json' \
  --data-binary @examples/field-g.json \
  http://localhost:8080/v1/run \
  --output result.h5
```

When authentication is configured, add
`-H "Authorization: Bearer $BELLHOP_AUTH_TOKEN"`.

## Errors

Errors use `application/json`:

```json
{
  "error": {
    "code": "invalid_case",
    "message": "request body describes an invalid case",
    "diagnostics": [
      {
        "severity": "error",
        "code": "BH0201",
        "message": "at least one launch angle is required",
        "field": "trace.launch_angles_degrees",
        "location": { "path": "request.json", "line": 1, "column": 1 }
      }
    ]
  }
}
```

| Status | Meaning |
| --- | --- |
| `400` | malformed JSON, wrong object shape, or unknown field/enum |
| `401` | configured Bearer token is missing or incorrect |
| `413` | request exceeds `BELLHOP_MAX_BODY_BYTES` |
| `415` | `Content-Type` is not JSON |
| `422` | semantic validation, unsupported run option, or simulation failure |
| `429` | a server-controlled simulation resource limit was exceeded |
| `500` | worker or HDF5 output failure |
| `504` | queueing plus execution exceeded the request timeout |

## Configuration

Every option is available as a command-line flag and environment variable.
Run `bellhop-server --help` for flag names.

| Environment variable | Default | Purpose |
| --- | ---: | --- |
| `BELLHOP_LISTEN` | `0.0.0.0:8080` | listener address |
| `BELLHOP_WORKERS` | logical CPU count | maximum admitted blocking requests |
| `BELLHOP_REQUEST_TIMEOUT_SECONDS` | `120` | queue plus execution timeout |
| `BELLHOP_MAX_BODY_BYTES` | `16777216` | request-body limit |
| `BELLHOP_AUTH_TOKEN` | unset | optional static Bearer token |
| `BELLHOP_MAX_RAYS` | `1000000` | total launch-ray limit |
| `BELLHOP_MAX_STEPS_PER_RAY` | `100000` | per-ray integration limit |
| `BELLHOP_MAX_TOTAL_POINTS` | `20000000` | ray/eigenray point limit |
| `BELLHOP_MAX_FIELD_CELLS` | `20000000` | pressure-field cell limit |

Logging uses `tracing` and the standard `RUST_LOG` filter; its fallback is
`info`.

CPU work runs in Tokio's blocking pool behind a Tower concurrency limit. A
request timeout returns `504` immediately, but Rust blocking work cannot be
forcefully cancelled: that simulation continues in the background until it
finishes and retains its worker slot. Size, concurrency, and solver limits
remain active. Cooperative solver cancellation
is intentionally deferred to a later version.

The optional token is deployment-wide, not per-user authorization. It is
compared on every `/v1/*` request; health and documentation remain public.
