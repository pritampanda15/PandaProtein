import React, { useEffect, useState } from 'react';
import { AlertCircle, CheckCircle2, Loader2 } from 'lucide-react';
import { api } from './api/client';
import { ClashesPanel } from './components/ClashesPanel';
import { CommandPalette } from './components/CommandPalette';
import { CootWorkbench } from './components/CootWorkbench';
import { ExportPanel } from './components/ExportPanel';
import { ImportModal } from './components/ImportModal';
import { LigandPanel } from './components/LigandPanel';
import { OverviewPanel } from './components/OverviewPanel';
import { PdbFixerPanel } from './components/PdbFixerPanel';
import { PropkaPanel } from './components/PropkaPanel';
import { RenameChainsPanel } from './components/RenameChainsPanel';
import { SequenceStrip } from './components/SequenceStrip';
import { SidebarNav } from './components/SidebarNav';
import { TopBar } from './components/TopBar';
import { ValidationPanel } from './components/ValidationPanel';
import { Viewer3D, type ChainHighlight } from './components/Viewer3D';
import type { NavTab, SelectedResidue, StructureSessionState } from './types/pdbprep';

interface Toast {
  id: number;
  type: 'success' | 'error';
  text: string;
}

export function App() {
  const [darkMode, setDarkMode] = useState(false);
  const [session, setSession] = useState<StructureSessionState | null>(null);
  const [activeTab, setActiveTab] = useState<NavTab>('overview');
  const [selectedResidue, setSelectedResidue] = useState<SelectedResidue | null>(null);
  const [focusTrigger, setFocusTrigger] = useState(0);
  const [highlightedChains, setHighlightedChains] = useState<ChainHighlight | null>(null);
  const [chainFocusTrigger, setChainFocusTrigger] = useState(0);
  const [viewerExpanded, setViewerExpanded] = useState(false);
  const [importOpen, setImportOpen] = useState(false);
  const [paletteOpen, setPaletteOpen] = useState(false);
  const [busy, setBusy] = useState<string | null>('Loading 6LU7 (SARS-CoV-2 Mpro)…');
  const [toast, setToast] = useState<Toast | null>(null);

  const pushToast = (text: string, type: 'success' | 'error' = 'success') => {
    const id = Date.now();
    setToast({ id, type, text });
    setTimeout(() => {
      setToast((prev) => (prev?.id === id ? null : prev));
    }, 3600);
  };

  // Sync Tailwind dark class on <html>
  useEffect(() => {
    document.documentElement.classList.toggle('dark', darkMode);
  }, [darkMode]);

  // Global macOS keyboard shortcuts: ⌘K Spotlight, ⌘Z Undo, ⌘⇧Z Redo
  useEffect(() => {
    const onKey = (e: KeyboardEvent) => {
      const mod = e.metaKey || e.ctrlKey;
      if (!mod) return;
      if (e.key.toLowerCase() === 'k') {
        e.preventDefault();
        setPaletteOpen((v) => !v);
      } else if (e.key.toLowerCase() === 'z' && !e.shiftKey) {
        const tag = (e.target as HTMLElement)?.tagName;
        if (tag === 'INPUT' || tag === 'TEXTAREA') return;
        if (session?.can_undo && !busy) {
          e.preventDefault();
          handleUndo();
        }
      } else if ((e.key.toLowerCase() === 'z' && e.shiftKey) || e.key.toLowerCase() === 'y') {
        const tag = (e.target as HTMLElement)?.tagName;
        if (tag === 'INPUT' || tag === 'TEXTAREA') return;
        if (session?.can_redo && !busy) {
          e.preventDefault();
          handleRedo();
        }
      }
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [session?.session_id, session?.can_undo, session?.can_redo, busy]);

  // Auto-load 6LU7 (SARS-CoV-2 Mpro + N3 inhibitor) as default structure on initial launch
  useEffect(() => {
    api
      .fetchRCSB('6LU7')
      .then((st) => {
        setSession(st);
        setBusy(null);
      })
      .catch((e) => {
        setBusy(null);
        pushToast(e?.message || 'Backend API unreachable — start FastAPI server on port 8000', 'error');
      });
  }, []);

  const handleSelectResidue = (res: SelectedResidue, focus3D = false) => {
    setSelectedResidue(res);
    if (focus3D) {
      setHighlightedChains(null);
      setFocusTrigger((n) => n + 1);
    }
  };

  const handleHighlightChains = (hl: ChainHighlight | null, zoom = false) => {
    setHighlightedChains(hl);
    if (zoom) {
      setChainFocusTrigger((n) => n + 1);
    }
  };

  const handleStateUpdate = (next: StructureSessionState, message: string) => {
    setSession(next);
    pushToast(message, 'success');
  };

  const handleUpload = async (file: File) => {
    setBusy(`Parsing ${file.name}…`);
    try {
      const st = await api.uploadFile(file);
      setSession(st);
      setSelectedResidue(null);
      setHighlightedChains(null);
      setFocusTrigger(0);
      setImportOpen(false);
      setActiveTab('overview');
      pushToast(`Loaded ${file.name} (${st.metadata.counts.protein_residues} residues)`);
    } catch (e: any) {
      pushToast(e?.message || 'Upload failed', 'error');
    } finally {
      setBusy(null);
    }
  };

  const handleFetchRCSB = async (pdbId: string) => {
    setBusy(`Downloading ${pdbId} from RCSB PDB…`);
    try {
      const st = await api.fetchRCSB(pdbId);
      setSession(st);
      setSelectedResidue(null);
      setHighlightedChains(null);
      setFocusTrigger(0);
      setImportOpen(false);
      setActiveTab('overview');
      pushToast(`Fetched RCSB ${pdbId} (${st.metadata.counts.protein_residues} residues)`);
    } catch (e: any) {
      pushToast(e?.message || `Failed to fetch ${pdbId}`, 'error');
    } finally {
      setBusy(null);
    }
  };

  const handleFetchAlphaFold = async (uniprotId: string) => {
    setBusy(`Fetching AlphaFold model ${uniprotId}…`);
    try {
      const st = await api.fetchAlphaFold(uniprotId);
      setSession(st);
      setSelectedResidue(null);
      setHighlightedChains(null);
      setFocusTrigger(0);
      setImportOpen(false);
      setActiveTab('overview');
      pushToast(`Fetched AlphaFold ${uniprotId}`);
    } catch (e: any) {
      pushToast(e?.message || `Failed to fetch AlphaFold ${uniprotId}`, 'error');
    } finally {
      setBusy(null);
    }
  };

  const handleLoadSample = async () => {
    setBusy('Loading Metalloenzyme Showcase…');
    try {
      const st = await api.loadSample('catalytic_showcase');
      setSession(st);
      setSelectedResidue(null);
      setHighlightedChains(null);
      setFocusTrigger(0);
      setImportOpen(false);
      setActiveTab('overview');
      pushToast('Loaded Metalloenzyme Diagnostic Showcase');
    } catch (e: any) {
      pushToast(e?.message || 'Failed to load sample', 'error');
    } finally {
      setBusy(null);
    }
  };

  const handleUndo = async () => {
    if (!session) return;
    setBusy('Undoing last edit…');
    try {
      const r = await api.undo(session.session_id);
      setSession(r.state);
      pushToast(r.undone ? `Undone: ${r.undone}` : 'Nothing to undo');
    } catch (e: any) {
      pushToast(e?.message || 'Undo failed', 'error');
    } finally {
      setBusy(null);
    }
  };

  const handleRedo = async () => {
    if (!session) return;
    setBusy('Redoing edit…');
    try {
      const r = await api.redo(session.session_id);
      setSession(r.state);
      pushToast(r.redone ? `Redone: ${r.redone}` : 'Nothing to redo');
    } catch (e: any) {
      pushToast(e?.message || 'Redo failed', 'error');
    } finally {
      setBusy(null);
    }
  };

  const handleDeepPrepPreset = async (
    preset: 'md_production' | 'docking_receptor' | 'crystallographic_polish' = 'md_production'
  ) => {
    if (!session) return;
    setBusy('Running Deep Prep (PDBFixer + COOT + PROPKA)…');
    try {
      const r = await api.deepPrepPreset(session.session_id, preset, session.propka.ph);
      setSession(r.state);
      pushToast(`Deep Prep complete — ${r.actions.length} automated repair steps applied`);
    } catch (e: any) {
      pushToast(e?.message || 'Deep Prep failed', 'error');
    } finally {
      setBusy(null);
    }
  };

  return (
    <div className="h-screen w-screen flex flex-col bg-[#F5F5F7] dark:bg-[#0D0D10] text-[#1D1D1F] dark:text-[#F5F5F7] overflow-hidden">
      <TopBar
        session={session}
        darkMode={darkMode}
        busy={busy}
        onToggleTheme={() => setDarkMode((d) => !d)}
        onOpenImport={() => setImportOpen(true)}
        onOpenCommandPalette={() => setPaletteOpen(true)}
        onUndo={handleUndo}
        onRedo={handleRedo}
        onQuickDeepPrep={() => handleDeepPrepPreset('md_production')}
        onGoExport={() => setActiveTab('export')}
      />

      <div className="flex-1 flex min-h-0 overflow-hidden">
        <SidebarNav activeTab={activeTab} onSelect={setActiveTab} session={session} />

        {/* Main Workspace Split View */}
        {!session ? (
          <div className="flex-1 flex flex-col items-center justify-center gap-3 p-8">
            <Loader2 className="w-7 h-7 animate-spin text-[#007AFF]" />
            <div className="text-[14px] font-medium text-[#6E6E73]">
              {busy || 'Loading structure workspace…'}
            </div>
            <button onClick={() => setImportOpen(true)} className="apple-btn-primary mt-2">
              Load Structure Manually
            </button>
          </div>
        ) : (
          <main className="flex-1 flex min-w-0 min-h-0 p-3 gap-3 overflow-hidden">
            {/* Left/Center Inspector Panel */}
            {!viewerExpanded && (
              <div className="w-full lg:w-[56%] xl:w-[54%] h-full overflow-y-auto pr-1">
                {activeTab === 'overview' && (
                  <OverviewPanel
                    session={session}
                    busy={!!busy}
                    onRunPreset={handleDeepPrepPreset}
                    onNavigate={setActiveTab}
                  />
                )}
                {activeTab === 'coot' && (
                  <CootWorkbench
                    session={session}
                    selectedResidue={selectedResidue}
                    busy={!!busy}
                    onSelectResidue={handleSelectResidue}
                    onStateUpdate={handleStateUpdate}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'chains' && (
                  <RenameChainsPanel
                    session={session}
                    busy={!!busy}
                    highlightedChains={highlightedChains}
                    onHighlightChains={handleHighlightChains}
                    onStateUpdate={handleStateUpdate}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'pdbfixer' && (
                  <PdbFixerPanel
                    session={session}
                    busy={!!busy}
                    onHighlightChains={handleHighlightChains}
                    onSelectResidue={handleSelectResidue}
                    onStateUpdate={handleStateUpdate}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'propka' && (
                  <PropkaPanel
                    session={session}
                    busy={!!busy}
                    onSelectResidue={handleSelectResidue}
                    onStateUpdate={handleStateUpdate}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'validation' && (
                  <ValidationPanel
                    session={session}
                    busy={!!busy}
                    onSelectResidue={handleSelectResidue}
                    onStateUpdate={handleStateUpdate}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'clashes' && (
                  <ClashesPanel
                    session={session}
                    busy={!!busy}
                    onSelectResidue={handleSelectResidue}
                    onStateUpdate={handleStateUpdate}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'ligands' && (
                  <LigandPanel
                    session={session}
                    busy={!!busy}
                    onSelectResidue={handleSelectResidue}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
                {activeTab === 'export' && (
                  <ExportPanel
                    session={session}
                    busy={!!busy}
                    onSetBusy={setBusy}
                    onError={(m) => pushToast(m, 'error')}
                  />
                )}
              </div>
            )}

            {/* Right Persistent Interactive Sequence Strip + 3D Molecular Viewport */}
            <div className={`${viewerExpanded ? 'w-full' : 'hidden lg:flex lg:w-[44%] xl:w-[46%]'} h-full min-h-0 flex flex-col gap-2.5`}>
              {session.sequence_chains && Object.keys(session.sequence_chains).length > 0 && (
                <SequenceStrip
                  sequenceChains={session.sequence_chains}
                  selectedResidue={selectedResidue}
                  onSelectResidue={(res, focus) => {
                    handleSelectResidue(res, focus);
                    if (activeTab === 'overview') setActiveTab('coot');
                  }}
                />
              )}
              <div className="flex-1 min-h-0">
                <Viewer3D
                  session={session}
                  darkMode={darkMode}
                  selectedResidue={selectedResidue}
                  onSelectResidue={(res) => {
                    handleSelectResidue(res, false);
                    if (activeTab === 'overview') {
                      setActiveTab('coot');
                    }
                  }}
                  focusTrigger={focusTrigger}
                  highlightedChains={highlightedChains}
                  chainFocusTrigger={chainFocusTrigger}
                  onClearChainHighlight={() => handleHighlightChains(null, true)}
                  expanded={viewerExpanded}
                  onToggleExpand={() => setViewerExpanded((v) => !v)}
                />
              </div>
            </div>
          </main>
        )}
      </div>

      {/* Import Modal */}
      <ImportModal
        open={importOpen}
        busy={!!busy}
        onClose={() => setImportOpen(false)}
        onUpload={handleUpload}
        onFetchRCSB={handleFetchRCSB}
        onFetchAlphaFold={handleFetchAlphaFold}
        onLoadSample={handleLoadSample}
      />

      {/* Apple Spotlight Command Palette (⌘K) */}
      <CommandPalette
        open={paletteOpen}
        session={session}
        onClose={() => setPaletteOpen(false)}
        onNavigate={setActiveTab}
        onSelectResidue={handleSelectResidue}
        onHighlightChains={handleHighlightChains}
        onQuickAction={(act) => {
          if (act === 'deep_prep') handleDeepPrepPreset('md_production');
          else if (act === 'coot_refine') handleDeepPrepPreset('crystallographic_polish');
          else if (act === 'open_import') setImportOpen(true);
          else if (act === 'toggle_theme') setDarkMode((d) => !d);
        }}
      />

      {/* Apple Floating Toast Notification */}
      {toast && (
        <div className="fixed bottom-5 left-1/2 -translate-x-1/2 z-50 pointer-events-none">
          <div
            className={`apple-glass px-4 py-2.5 rounded-2xl shadow-apple-lg flex items-center gap-2 text-[13px] font-medium ${
              toast.type === 'error'
                ? 'border-[#FF3B30]/40 text-[#D70015] dark:text-[#FF453A]'
                : 'border-[#34C759]/35 text-[#1D1D1F] dark:text-white'
            }`}
          >
            {toast.type === 'error' ? (
              <AlertCircle className="w-4 h-4 text-[#FF3B30] shrink-0" />
            ) : (
              <CheckCircle2 className="w-4 h-4 text-[#34C759] shrink-0" />
            )}
            <span>{toast.text}</span>
          </div>
        </div>
      )}
    </div>
  );
}

export default App;
