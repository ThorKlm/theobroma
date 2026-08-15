-- Supplement to apg_clade_ref: fills the 2,944-compound unmapped tail.
-- Group A: genuine angiosperm families omitted from the first pass (correct APG IV clade).
-- Group B: NON-angiosperm families that are misclassified as magnoliopsida/liliopsida in
--          the source taxonomy; labeled by true higher group so apg_clade is honest and
--          the misclassification is visible (these are a data-quality finding).
-- Idempotent via ON CONFLICT.

INSERT INTO apg_clade_ref (family, clade) VALUES
-- ===== Group A: genuine angiosperms (eudicots) =====
('calophyllaceae','eudicot'),('gelsemiaceae','eudicot'),('lardizabalaceae','eudicot'),
('ancistrocladaceae','eudicot'),('balanophoraceae','eudicot'),('hydrangeaceae','eudicot'),
('melastomataceae','eudicot'),('aphloiaceae','eudicot'),('olacaceae','eudicot'),
('staphyleaceae','eudicot'),('dioncophyllaceae','eudicot'),('nyssaceae','eudicot'),
('loasaceae','eudicot'),('peraceae','eudicot'),('petiveriaceae','eudicot'),
('francoaceae','eudicot'),('stachyuraceae','eudicot'),('cynomoriaceae','eudicot'),
('vochysiaceae','eudicot'),('crossosomataceae','eudicot'),('loasaceae','eudicot'),
('aextoxicaceae','eudicot'),('berberidopsidaceae','eudicot'),('haptanthaceae','eudicot'),
('petenaeaceae','eudicot'),('nothofagaceae','eudicot'),('ticodendraceae','eudicot'),
('geissolomataceae','eudicot'),('strasburgeriaceae','eudicot'),('vivianiaceae','eudicot'),
('ledocarpaceae','eudicot'),('melianthaceae','eudicot'),('francoaceae','eudicot'),
('greyiaceae','eudicot'),('rhipogonaceae','monocot'),('ripogonaceae','monocot'),
('petermanniaceae','monocot'),('philesiaceae','monocot'),('behniaceae','monocot'),
('herreriaceae','monocot'),('aphyllanthaceae','monocot'),('acoraceae','monocot'),
('maundiaceae','monocot'),('posidoniaceae','monocot'),('ruppiaceae','monocot'),
('cymodoceaceae','monocot'),('zannichelliaceae','monocot'),('najadaceae','monocot'),
('limnocharitaceae','monocot'),
-- ===== Group A: magnoliids omitted =====
('eupteleaceae','eudicot'),  -- Eupteleaceae is Ranunculales (eudicot), not magnoliid
-- ===== Group B: NON-angiosperms misclassified under magnoliopsida/liliopsida =====
-- Gymnosperms (conifers, cycads, gnetophytes)
('cupressaceae','gymnosperm'),('pinaceae','gymnosperm'),('araucariaceae','gymnosperm'),
('taxaceae','gymnosperm'),('podocarpaceae','gymnosperm'),('ephedraceae','gymnosperm'),
('cephalotaxaceae','gymnosperm'),('sciadopityaceae','gymnosperm'),('cycadaceae','gymnosperm'),
('zamiaceae','gymnosperm'),('ginkgoaceae','gymnosperm'),('gnetaceae','gymnosperm'),
('welwitschiaceae','gymnosperm'),('taxodiaceae','gymnosperm'),
-- Ferns and lycophytes (monilophytes / lycophytes)
('aspleniaceae','fern-lycophyte'),('pteridaceae','fern-lycophyte'),
('polypodiaceae','fern-lycophyte'),('lycopodiaceae','fern-lycophyte'),
('equisetaceae','fern-lycophyte'),('dennstaedtiaceae','fern-lycophyte'),
('salviniaceae','fern-lycophyte'),('osmundaceae','fern-lycophyte'),
('lygodiaceae','fern-lycophyte'),('schizaeaceae','fern-lycophyte'),
('dipteridaceae','fern-lycophyte'),('selaginellaceae','fern-lycophyte'),
('dryopteridaceae','fern-lycophyte'),('blechnaceae','fern-lycophyte'),
('thelypteridaceae','fern-lycophyte'),('athyriaceae','fern-lycophyte'),
('cyatheaceae','fern-lycophyte'),('gleicheniaceae','fern-lycophyte'),
('marattiaceae','fern-lycophyte'),('ophioglossaceae','fern-lycophyte'),
('psilotaceae','fern-lycophyte'),('isoetaceae','fern-lycophyte'),
('marsileaceae','fern-lycophyte'),('hymenophyllaceae','fern-lycophyte'),
('davalliaceae','fern-lycophyte'),('nephrolepidaceae','fern-lycophyte'),
('woodsiaceae','fern-lycophyte'),('onocleaceae','fern-lycophyte'),
-- Bryophytes (liverworts, hornworts, mosses)
('anthocerotaceae','bryophyte'),('aneuraceae','bryophyte'),('herbertaceae','bryophyte'),
('monocleaceae','bryophyte'),('lepidoziaceae','bryophyte'),('lophocoleaceae','bryophyte'),
('polytrichaceae','bryophyte'),('marchantiaceae','bryophyte'),('ricciaceae','bryophyte'),
('jungermanniaceae','bryophyte'),('scapaniaceae','bryophyte'),('porellaceae','bryophyte'),
('frullaniaceae','bryophyte'),('radulaceae','bryophyte'),('plagiochilaceae','bryophyte'),
('sphagnaceae','bryophyte'),('polytrichaceae','bryophyte'),('bryaceae','bryophyte'),
('dicranaceae','bryophyte'),('hypnaceae','bryophyte'),('mniaceae','bryophyte'),
('conocephalaceae','bryophyte'),('lunulariaceae','bryophyte'),('targioniaceae','bryophyte'),
-- Algae (misclassified as land plants)
('rhodomelaceae','algae'),('caulerpaceae','algae'),('bangiaceae','algae'),
('corallinaceae','algae'),('gracilariaceae','algae'),('sargassaceae','algae'),
('dictyotaceae','algae'),('ulvaceae','algae'),('codiaceae','algae'),
('gigartinaceae','algae'),('ceramiaceae','algae'),('cladophoraceae','algae')
ON CONFLICT (family) DO NOTHING;

-- corynocarpaceae is Cucurbitales (eudicot); trimeniaceae is basal-angiosperm;
-- gomortegaceae is magnoliid; tovariaceae, tropaeolaceae, griseliniaceae, calyceraceae,
-- montiaceae, vahliaceae, lentibulariaceae, irvingiaceae, biebersteiniaceae are eudicots
INSERT INTO apg_clade_ref (family, clade) VALUES
('corynocarpaceae','eudicot'),('trimeniaceae','basal-angiosperm'),
('gomortegaceae','magnoliid'),('tovariaceae','eudicot'),('tropaeolaceae','eudicot'),
('griseliniaceae','eudicot'),('calyceraceae','eudicot'),('montiaceae','eudicot'),
('vahliaceae','eudicot'),('lentibulariaceae','eudicot'),('irvingiaceae','eudicot'),
('biebersteiniaceae','eudicot'),('marantaceae','monocot')
ON CONFLICT (family) DO NOTHING;
