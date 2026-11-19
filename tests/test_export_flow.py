"""End-to-end tests of the Export tab: "Download all (.zip)" with pending GSEA collections
and the "Run pending GSEA now" button.

They start the real Shiny app on a local port, speak Shiny's websocket protocol like a
browser (load the example data, set inputs, click buttons) and download the zip over
HTTP, so they cover the part that used to break: work done inside the streaming
download after Shiny had already sent "200 OK" (0-byte / 33-byte zips).

GSEA is kept fast by limiting the run to two collections with few permutations.
"""

import asyncio
import io
import json
import re
import socket
import threading
import time
import urllib.request
import zipfile
from html.parser import HTMLParser

import pytest
import uvicorn
import websockets

import ezplot.app as ezapp
from ezplot import genesets
from ezplot.project import GseaSettings, Project

KEYS = ["hallmarks", "go_cc"]


class FastProject(Project):
    """Example-sized GSEA: two collections, few permutations."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.gsea = GseaSettings(nperm=200)

    def collections(self):
        return {k: v for k, v in super().collections().items() if k in KEYS}


class FailingProject(FastProject):
    """GO CC GSEA raises when run, and the heatmap figure raises."""

    def run_gsea(self, key, progress=None):
        if key == "go_cc" and not self.gsea_done(key):
            raise RuntimeError("boom in GSEA")
        return super().run_gsea(key, progress)

    def heatmap_fig(self, p=None):
        raise ValueError("boom in heatmap")


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture
def server(monkeypatch, request):
    monkeypatch.setattr(ezapp, "Project", request.param)
    port = _free_port()
    srv = uvicorn.Server(uvicorn.Config(ezapp.app, host="127.0.0.1", port=port,
                                        log_level="warning", ws_ping_interval=None))
    th = threading.Thread(target=srv.run, daemon=True)
    th.start()
    for _ in range(200):
        if srv.started:
            break
        time.sleep(0.05)
    yield port
    srv.should_exit = True
    th.join(10)


class _Inputs(HTMLParser):
    """Default values of every input on the page, as the browser would send them."""

    def __init__(self):
        super().__init__()
        self.values, self._select = {}, None

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        typ, id_ = a.get("type"), a.get("id")
        if tag == "input" and id_ and typ in ("number", "text", "color"):
            v = a.get("value") or ""
            self.values[id_] = (float(v) if "." in v else int(v)) if typ == "number" and v else v
        elif tag == "input" and typ == "checkbox" and id_:
            self.values[id_] = "checked" in a
        elif tag == "input" and typ == "radio" and "checked" in a:
            self.values[a["name"]] = a["value"]
        elif tag == "input" and typ == "checkbox" and a.get("name"):
            group = self.values.setdefault(a["name"], [])
            if "checked" in a:
                group.append(a["value"])
        elif tag == "textarea" and id_:
            self.values[id_] = ""
        elif tag == "select" and id_:
            self._select = id_
            self.values[id_] = None
        elif tag == "option" and self._select:
            if self.values[self._select] is None or "selected" in a:
                self.values[self._select] = a.get("value")

    def handle_endtag(self, tag):
        if tag == "select":
            self._select = None


class Client:
    def __init__(self, port):
        self.port, self.ws, self.outputs, self.notifications = port, None, {}, []

    async def __aenter__(self):
        html = urllib.request.urlopen(f"http://127.0.0.1:{self.port}/").read().decode()
        p = _Inputs()
        p.feed(html)
        inputs = {**p.values, "tabs": "Export", "example:shiny.action": 0,
                  "gsea_apply:shiny.action": 0, "group_apply:shiny.action": 0,
                  "x_run_pending:shiny.action": 0, "g_pathways": None, "files": None,
                  ".clientdata_url_search": "", ".clientdata_url_hash_initial": "",
                  # outputs are suspended until the browser reports them visible
                  ".clientdata_output_x_status_hidden": False}
        self.ws = await websockets.connect(f"ws://127.0.0.1:{self.port}/websocket/",
                                           max_size=None, ping_interval=None)
        self.session_id = json.loads(await self.ws.recv())["config"]["sessionId"]
        await self.ws.send(json.dumps({"method": "init", "data": inputs}))
        self.reader = asyncio.create_task(self._read())
        return self

    async def __aexit__(self, *exc):
        self.reader.cancel()
        await self.ws.close()

    async def _read(self):
        async for raw in self.ws:
            msg = json.loads(raw)
            self.outputs.update(msg.get("values") or {})
            note = msg.get("notification")
            if note and note.get("type") == "show":
                self.notifications.append(note["message"])

    async def update(self, **values):
        await self.ws.send(json.dumps({"method": "update", "data": values}))

    async def wait_for(self, cond, timeout=120):
        end = time.time() + timeout
        while time.time() < end:
            if cond():
                return
            await asyncio.sleep(0.1)
        raise AssertionError("timed out waiting for the app")

    def status(self):
        out = self.outputs.get("x_status")
        return re.sub(r"<[^>]+>", "", out["html"] if isinstance(out, dict) else out or "")

    async def download(self, name):
        url = f"http://127.0.0.1:{self.port}/session/{self.session_id}/download/{name}?w="

        def get():
            with urllib.request.urlopen(url, timeout=300) as r:
                return r.status, r.read()

        return await asyncio.get_running_loop().run_in_executor(None, get)

    async def load_example(self):
        await self.update(**{"example:shiny.action": 1})
        await self.wait_for(lambda: "Not run with the current settings" in self.status())


def _pending(status):
    return status.split("Not run with the current settings:")[1].strip()


def _zip(data):
    assert len(data) > 1000, f"broken zip of {len(data)} bytes: {data[:80]!r}"
    z = zipfile.ZipFile(io.BytesIO(data))
    assert z.testzip() is None
    return z


def run(coro):
    return asyncio.run(asyncio.wait_for(coro, 600))


@pytest.mark.parametrize("server", [FastProject], indirect=True)
def test_download_all_runs_pending_gsea(server):
    async def go():
        async with Client(server) as c:
            await c.load_example()
            assert _pending(c.status()) != "none"
            await c.update(x_runall=True)
            await asyncio.sleep(0.3)
            status, data = await c.download("dl_all")
            assert status == 200
            names = _zip(data).namelist()
            for key in KEYS:
                assert f"tables/gsea_{key}.tsv" in names
                assert f"figures/gsea_{key}.png" in names
            assert "figures/volcano.png" in names and "figures/heatmap.png" in names
            assert "errors.txt" not in names
            # x_status was refreshed by the download itself
            await c.wait_for(lambda: _pending(c.status()) == "none", 30)
    run(go())


@pytest.mark.parametrize("server", [FailingProject], indirect=True)
def test_download_all_with_failures_still_gives_valid_zip(server):
    async def go():
        async with Client(server) as c:
            await c.load_example()
            await c.update(x_runall=True)
            await asyncio.sleep(0.3)
            status, data = await c.download("dl_all")
            assert status == 200
            z = _zip(data)
            names = z.namelist()
            assert "errors.txt" in names
            err = z.read("errors.txt").decode()
            assert "boom in GSEA" in err and "boom in heatmap" in err
            assert "figures/volcano.png" in names           # the rest is still there
            assert "tables/gsea_hallmarks.tsv" in names
            assert "tables/gsea_go_cc.tsv" not in names
            assert "ezplot_parameters.json" in names
    run(go())


@pytest.mark.parametrize("server", [FastProject], indirect=True)
def test_run_pending_button_updates_status(server):
    async def go():
        async with Client(server) as c:
            await c.load_example()
            assert _pending(c.status()) != "none"
            await c.update(**{"x_run_pending:shiny.action": 1})
            await c.wait_for(lambda: _pending(c.status()) == "none", 120)
            for title in (genesets.COLLECTIONS[k][0].replace("GSEA ", "") for k in KEYS):
                assert title in c.status().split("Not run")[0]
            # with everything run, the zip has all collections without x_runall
            status, data = await c.download("dl_all")
            names = _zip(data).namelist()
            assert all(f"tables/gsea_{k}.tsv" in names for k in KEYS)
    run(go())


@pytest.mark.parametrize("server", [FailingProject], indirect=True)
def test_run_pending_button_reports_errors(server):
    async def go():
        async with Client(server) as c:
            await c.load_example()
            await c.update(**{"x_run_pending:shiny.action": 1})
            await c.wait_for(lambda: any("boom in GSEA" in json.dumps(n)
                                         for n in c.notifications), 120)
            # hallmarks still ran and the session is alive
            await c.wait_for(lambda: "Hallmarks" in c.status().split("Not run")[0], 30)
            assert "GO Cellular Component" in _pending(c.status())
    run(go())
