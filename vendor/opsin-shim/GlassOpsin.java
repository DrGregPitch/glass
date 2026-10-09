package uk.ac.cam.ch.wwmm.opsin;
// Glass sidecar around OPSIN's own parse pipeline.  For a name it reports the raw token tree
// (before OPSIN rewrites locants into attributes, so character offsets can be recovered) and,
// for every token that owns a fragment, the atom ids after group generation ("pre") and after
// suffix application ("post").  Name regions then come from OPSIN's grammar, not a vocabulary.
// Radicals are allowed, so substituent names (methyl, phenyl) parse too.
// Protocol: one name per stdin line; one JSON object per stdout line.
import java.io.*;
import java.lang.reflect.Field;
import java.util.*;

public class GlassOpsin {
    public static void main(String[] args) throws Exception {
        NameToStructure nts = NameToStructure.getInstance();
        Field pf = NameToStructure.class.getDeclaredField("parser"); pf.setAccessible(true);
        Field sf = NameToStructure.class.getDeclaredField("suffixRules"); sf.setAccessible(true);
        Parser parser = (Parser) pf.get(nts);
        SuffixRules suffixRules = (SuffixRules) sf.get(nts);
        BufferedReader in = new BufferedReader(new InputStreamReader(System.in, "UTF-8"));
        PrintStream out = new PrintStream(new FileOutputStream(FileDescriptor.out), true, "UTF-8");
        out.println("{\"ready\":true,\"opsin\":" + q(NameToStructure.getVersion()) + "}");
        String line;
        while ((line = in.readLine()) != null) {
            String name = line.trim();
            if (!name.isEmpty()) out.println(handle(name, parser, suffixRules));
        }
    }

    static String handle(String name, Parser parser, SuffixRules suffixRules) {
        String err = "";
        try {
            NameToStructureConfig config = NameToStructureConfig.getDefaultConfigInstance().clone();
            config.setAllowRadicals(true);
            String pre = PreProcessor.preProcess(name);
            List<Element> parses = parser.parse(config, pre);
            Collections.sort(parses, new SortParses());
            for (Element parse : parses) {
                try {
                    // raw tree first: every token in name order, keyed by identity
                    List<Element> order = new ArrayList<Element>();
                    Map<Element, Integer> idx = new IdentityHashMap<Element, Integer>();
                    List<Integer> parents = new ArrayList<Integer>();
                    collect(parse, -1, order, idx, parents);
                    BuildState state = new BuildState(config);
                    new ComponentGenerator(state).processParse(parse);
                    new ComponentProcessor(state, new SuffixApplier(state, suffixRules)).processParse(parse);
                    Map<Integer, String> postAtoms = atomsOf(order, idx);
                    // multiplied substituents (dimethyl) are cloned by OPSIN: report the clones' atoms too
                    Fragment frag = new StructureBuilder(state).buildFragment(parse);
                    StringBuilder extras = new StringBuilder();
                    List<Element> after = new ArrayList<Element>(); collect(parse, -1, after, new IdentityHashMap<Element, Integer>(), new ArrayList<Integer>());
                    for (Element el : after) {
                        if (idx.containsKey(el) || !(el instanceof TokenEl)) continue;
                        Fragment f; try { f = el.getFrag(); } catch (Exception e) { continue; }
                        if (f == null) continue;
                        StringBuilder ids = new StringBuilder();
                        for (Atom at : f.getAtomList()) { if (ids.length() > 0) ids.append(','); ids.append(at.getID()); }
                        if (extras.length() > 0) extras.append(',');
                        extras.append("{\"el\":").append(q(el.getName())).append(",\"v\":").append(q(el.getValue())).append(",\"atoms\":[").append(ids).append("]}");
                    }
                    OpsinResult r = new OpsinResult(frag, OpsinResult.OPSIN_RESULT_STATUS.SUCCESS, "", name);
                    StringBuilder toks = new StringBuilder();
                    for (int i = 0; i < order.size(); i++) {
                        Element el = order.get(i);
                        boolean leaf = el instanceof TokenEl;
                        if (toks.length() > 0) toks.append(',');
                        toks.append("{\"i\":").append(i).append(",\"p\":").append(parents.get(i))
                            .append(",\"el\":").append(q(el.getName())).append(",\"v\":").append(q(leaf ? el.getValue() : null))
                            .append(",\"g\":").append(q(leaf ? el.getAttributeValue("value") : null))
                            .append(",\"atoms\":[").append(postAtoms.containsKey(i) ? postAtoms.get(i) : "").append("]}");
                    }
                    return "{\"name\":" + q(name) + ",\"pre\":" + q(pre) + ",\"smiles\":" + q(r.getSmiles()) + ",\"cml\":" + q(r.getCml())
                         + ",\"tokens\":[" + toks + "],\"extras\":[" + extras + "]}";
                } catch (Exception e) { err = String.valueOf(e.getMessage()); }
            }
        } catch (Exception e) { err = String.valueOf(e.getMessage()); }
        return "{\"name\":" + q(name) + ",\"error\":" + q(err.isEmpty() ? "no parse" : err) + "}";
    }

    static void collect(Element el, int parent, List<Element> order, Map<Element, Integer> idx, List<Integer> parents) {
        int me = order.size(); idx.put(el, me); order.add(el); parents.add(parent);
        for (Element c : el.getChildElements()) collect(c, me, order, idx, parents);
    }

    static Map<Integer, String> atomsOf(List<Element> order, Map<Element, Integer> idx) {
        Map<Integer, String> m = new HashMap<Integer, String>();
        for (Element el : order) {
            if (!(el instanceof TokenEl)) continue;
            Fragment f;
            try { f = el.getFrag(); } catch (Exception e) { continue; }
            if (f == null) continue;
            StringBuilder ids = new StringBuilder();
            for (Atom at : f.getAtomList()) { if (ids.length() > 0) ids.append(','); ids.append(at.getID()); }
            m.put(idx.get(el), ids.toString());
        }
        return m;
    }

    static String q(String s) {
        if (s == null) return "null";
        StringBuilder b = new StringBuilder("\"");
        for (char c : s.toCharArray()) {
            if (c == '"' || c == '\\') b.append('\\').append(c);
            else if (c < 0x20) b.append(String.format("\\u%04x", (int) c));
            else b.append(c);
        }
        return b.append('"').toString();
    }
}
