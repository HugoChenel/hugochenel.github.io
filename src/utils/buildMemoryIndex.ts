export interface MemoryEntry {
  id: string;
  title: string;
  url: string;
  category: string;
  excerpt: string;
  tags: string[];
  date?: string;
  headings: string[];
  body: string;
}

type MarkdownModule = {
  frontmatter?: Record<string, unknown>;
  rawContent?: () => string;
};

// ─────────────────────────────────────────────────────────────────────────────
// Content cleaning — applied at build time so the browser index contains
// readable natural text only. Never called client-side.
// ─────────────────────────────────────────────────────────────────────────────
export function cleanContentForIndexing(raw: string): string {
  return (
    raw
      // Fenced code blocks (may be multiline)
      .replace(/```[\s\S]*?```/g, ' ')
      // Inline code
      .replace(/`[^`\n]{0,200}`/g, ' ')
      // HTML comments
      .replace(/<!--[\s\S]*?-->/g, ' ')
      // HTML/JSX tags (self-closing and paired, up to 600 chars to avoid catastrophic backtracking)
      .replace(/<[^>]{0,600}>/g, ' ')
      // Markdown image syntax
      .replace(/!\[[^\]]{0,200}\]\([^)]{0,500}\)/g, ' ')
      // Markdown links → keep link text
      .replace(/\[([^\]]{0,200})\]\([^)]{0,500}\)/g, '$1')
      // Markdown headings → keep heading text (strip # markers)
      .replace(/^#{1,6}\s+/gm, '')
      // import / export statements
      .replace(/^(import|export)\s.*$/gm, ' ')
      // Blockquotes
      .replace(/^>\s*/gm, '')
      // Horizontal rules
      .replace(/^[-*_]{3,}\s*$/gm, ' ')
      // Bold / italic / strikethrough
      .replace(/\*{1,3}([^*\n]{0,300})\*{1,3}/g, '$1')
      .replace(/_{1,3}([^_\n]{0,300})_{1,3}/g, '$1')
      .replace(/~~([^~\n]{0,300})~~/g, '$1')
      // Bare URLs
      .replace(/https?:\/\/\S+/g, ' ')
      // Table pipes
      .replace(/\|/g, ' ')
      // Normalise whitespace
      .replace(/\s+/g, ' ')
      .trim()
  );
}

export function buildBlogEntries(
  modules: Record<string, MarkdownModule>
): MemoryEntry[] {
  return Object.entries(modules).map(([path, mod]) => {
    const fm = (mod.frontmatter ?? {}) as Record<string, unknown>;
    const slug = (path.split('/').pop() ?? '').replace(/\.(md|mdx)$/, '');

    let rawBody = '';
    try {
      rawBody = typeof mod.rawContent === 'function' ? mod.rawContent() : '';
    } catch {
      // rawContent unavailable in this build mode
    }

    // Extract headings from raw body before cleaning
    const headings = Array.from(
      rawBody.matchAll(/^#{1,3} +(.+)$/gm),
      (m) => m[1].trim()
    );

    // Clean at build time — browser only ever sees readable text
    const cleanBody = cleanContentForIndexing(rawBody);

    const body = [
      String(fm.title ?? ''),
      String(fm.excerpt ?? ''),
      ((fm.tags as string[] | undefined) ?? []).join(' '),
      cleanBody.slice(0, 2500),
    ]
      .join(' ')
      .trim();

    return {
      id: `blog-${slug}`,
      title: String(fm.title ?? slug),
      url: `/blog/${slug}`,
      category: String((fm.tags as string[] | undefined)?.[0] ?? 'Blog'),
      excerpt: String(fm.excerpt ?? ''),
      tags: (fm.tags as string[] | undefined) ?? [],
      date: fm.date != null ? String(fm.date) : undefined,
      headings,
      body,
    };
  });
}

export const CURATED_PAGE_ENTRIES: MemoryEntry[] = [
  {
    id: 'page-home',
    title: 'Hugo Chenel',
    url: '/',
    category: 'About',
    excerpt: 'PhD researcher in computational oncology, hybrid athlete, and biohacker based in Toulouse.',
    tags: ['about', 'Hugo Chenel', 'science', 'biohacking', 'athleticism', 'Toulouse'],
    headings: ['Science', 'Biohacking', 'Athleticism'],
    body: 'Hugo Chenel is a PhD researcher in computational oncology at the CRCT in Toulouse, France. He combines three domains: rigorous scientific research on cancer, systematic biohacking for human optimisation, and high-performance hybrid athletics. He competes in Hyrox, CrossFit Open, and road races. His research applies machine learning to transcriptomics data to model cancer cell dynamics.',
  },
  {
    id: 'page-science',
    title: 'Science & Research',
    url: '/science',
    category: 'Science',
    excerpt: 'PhD research in computational oncology, ML models of cancer dynamics, single-cell transcriptomics.',
    tags: [
      'science', 'research', 'PhD', 'cancer', 'machine learning', 'bioinformatics',
      'transcriptomics', 'oncology', 'RNA-seq', 'single-cell', 'scRNA-seq', 'DESeq2',
      'computational biology', 'CRCT', 'Toulouse', 'tumor', 'gene expression',
    ],
    headings: ['Research focus', 'Publications', 'PhD Thesis', 'Tutorials', 'Conferences'],
    body: `Machine learning models help study cancer by identifying patterns in high-dimensional gene expression data that are too complex for manual analysis. Single-cell RNA-seq (scRNA-seq) measures the transcriptome of individual cells, revealing tumor heterogeneity that bulk sequencing misses. Dimensionality reduction methods like PCA and UMAP, followed by clustering algorithms, group cells into subtypes based on expression profiles. Trajectory inference reconstructs how cancer cells transition between states over time. Computational oncology combines these methods to discover therapeutic targets and understand treatment resistance.

PhD thesis: "Machine learning models of cancer dynamics based on transcriptomics data", supervised by Dr. Vera Pancaldi and Dr. Andrei Zinovyev at the CRCT (Centre de Recherches en Cancérologie de Toulouse), funded by Toulouse University, completion expected 2026.

Bulk RNA-seq analysis uses DESeq2 in R to identify differentially expressed genes between conditions. The pipeline covers quality control, alignment, count normalization, and downstream functional enrichment. Single-cell RNA-seq preprocessing with Scanpy (Python) or Seurat (R) includes filtering low-quality cells, normalization, feature selection, dimensionality reduction, clustering, and cell type annotation.`,
  },
  {
    id: 'page-biohacking',
    title: 'Biohacking',
    url: '/biohacking',
    category: 'Health',
    excerpt: 'Systematic health optimisation: sleep, supplements, biomarkers, and advanced protocols.',
    tags: [
      'biohacking', 'health', 'supplements', 'sleep', 'protocol', 'biomarkers',
      'HRV', 'testosterone', 'HBOT', 'cold', 'Whoop', 'creatine', 'omega-3',
      'magnesium', 'vitamin D', 'AG1', 'float tank', 'sensory deprivation',
      'recovery', 'longevity', 'optimisation',
    ],
    headings: [
      'Master the basics', 'First principles', 'Advanced protocol',
      'Supplement stack', 'Biomarkers', 'Genetics',
    ],
    body: `Biohacking is the systematic approach to optimising human biology through data, protocols, and evidence-based interventions. The foundation is mastering basics: consistent sleep timing, whole-food nutrition, daily movement, stress management, and social connection.

Supplement stack: AG1 (comprehensive greens), creatine monohydrate 5g daily (proven for strength and cognition), omega-3 fish oil (anti-inflammatory, cardiovascular), magnesium glycinate (sleep and muscle recovery), vitamin D3 with K2 (immune and bone health), zinc (testosterone and immune support). Creatine is one of the most evidence-backed supplements for athletic performance and cognitive function.

Biomarkers tracked: testosterone 850 ng/dL, free testosterone 26.2 pg/mL. HRV monitored daily via Whoop band to gauge recovery status and training readiness.

Advanced protocols: Hyperbaric oxygen therapy (HBOT) — breathing pure oxygen at 1.5–2.5 atm increases oxygen delivery to tissues, accelerates healing, reduces inflammation, and may enhance cognitive performance. Sensory deprivation float tanks eliminate external stimulation to induce deep relaxation and promote neural recovery. Cold exposure triggers hormetic adaptation: reduces inflammation, activates brown fat, and may boost testosterone.

Supplements support athletic recovery by reducing muscle damage, accelerating protein synthesis, and maintaining hormonal balance. Magnesium supports over 300 enzymatic processes and is frequently depleted in athletes. Omega-3 reduces exercise-induced inflammation and supports cardiovascular health.`,
  },
  {
    id: 'page-athleticism',
    title: 'Athleticism & Training',
    url: '/athleticism',
    category: 'Fitness',
    excerpt: 'Hybrid athlete competing in Hyrox, CrossFit, and road races. VO2max, PRs, weekly training.',
    tags: [
      'athleticism', 'Hyrox', 'CrossFit', 'training', 'fitness', 'strength',
      'VO2max', 'hybrid athlete', 'running', 'racing', 'Red Bull 400',
      'zone-2', 'cardio', 'aerobic', 'HIIT', 'deadlift', 'squat', 'bench',
      'endurance', 'performance',
    ],
    headings: [
      'Foundations', 'Physiological profile', 'Competition results',
      'Personal records', 'Weekly training schedule',
    ],
    body: `Athleticism here covers hybrid performance — combining endurance and strength across Hyrox, CrossFit Open, and road races. Hyrox is a global fitness race combining 8 km of running with 8 functional strength stations: SkiErg, sled push, sled pull, burpee broad jumps, rowing, farmers carry, sandbag lunges, and wall balls. Training for Hyrox requires both aerobic capacity and muscular endurance. Effective preparation includes zone-2 cardio for aerobic base, threshold runs for pace conditioning, and station-specific strength work. Best Hyrox result: 01:07:55, Toulouse, March 2026, Open Solo.

Zone-2 cardio is sustained aerobic exercise at roughly 60–70% of maximum heart rate — the intensity where you can hold a conversation. It is the most effective training zone for building mitochondrial density, improving fat oxidation, and developing a durable aerobic base.

VO2max is the maximum rate of oxygen consumption during exercise, expressed in mL/kg/min, and is the single best predictor of cardiovascular fitness. To improve VO2max: combine zone-2 base training with high-intensity intervals. Physiological profile: VO2max 62 mL/kg/min, resting heart rate ~42 bpm, body fat 8.5–9.5%.

Hybrid athlete training weekly structure: zone-2 cardio runs (3–4×/week), strength training focused on compound lifts (3×/week), one HIIT or Hyrox-specific session, and mobility work. Strength PRs: bench press 102.5 kg, back squat 110 kg, deadlift 160 kg, overhead press 75 kg. CrossFit Open 2026: Level 8, 13–23 percentile, Advanced category.`,
  },
  {
    id: 'page-outreach',
    title: 'Outreach',
    url: '/outreach',
    category: 'Outreach',
    excerpt: 'Media appearances and public communication around cancer research, sport, and health in France.',
    tags: ['outreach', 'media', 'science communication', 'cancer', 'public speaking', 'France'],
    headings: ['Apparitions Médias', 'Mission', 'Prévention', 'Recherche'],
    body: `Science outreach and public communication around cancer research, athletic performance, and health optimisation. Media appearances in France covering PhD research in computational oncology, cancer prevention, and the intersection of sport and science. Cancer prevention: 40% of cancers can be avoided through better information and lifestyle choices. The outreach mission covers prevention education, patient support, and research funding advocacy.`,
  },
  {
    id: 'page-skydiving',
    title: 'Skydiving',
    url: '/blog/Skydiving',
    category: 'Adventure',
    excerpt: 'Licensed skydiver — the adrenaline dimension of a performance-focused life.',
    tags: ['skydiving', 'adventure', 'adrenaline', 'sport'],
    headings: ['Skydiving'],
    body: `Skydiving is a licensed sport involving freefall from aircraft at altitude. It develops spatial awareness, stress inoculation, and comfort with controlled risk — transferable mental skills for high-performance sport and research environments.`,
  },
];
