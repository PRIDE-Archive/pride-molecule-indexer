package uk.ac.ebi.pride.archive.indexer.utility;

import org.ehcache.Cache;
import org.ehcache.CacheManager;
import org.ehcache.config.builders.CacheConfigurationBuilder;
import org.ehcache.config.builders.CacheManagerBuilder;
import org.ehcache.config.builders.ResourcePoolsBuilder;
import org.ehcache.config.units.EntryUnit;
import org.ehcache.config.units.MemoryUnit;

import java.io.File;
import java.io.Serializable;
import java.util.List;

public class AppCacheManager implements Serializable {

    public final static long serialVersionUID = 1L;
    private static AppCacheManager instance = null;
    private final static CacheManager cacheManage;

    private static final String INPUT_SPECTRA_CACHE = "InputSpectra";
    private static final String CLUSTERS_CACHE = "Clusters";
    private static final String PROTEIN_TO_PSMS_CACHE = "ProteinToPsms";

    static {
        ResourcePoolsBuilder spectraResourceBuilder = ResourcePoolsBuilder.newResourcePoolsBuilder()
                .heap(100000, EntryUnit.ENTRIES)
                .offheap(300, MemoryUnit.MB)
                .disk(5, MemoryUnit.GB);

        ResourcePoolsBuilder proteinsResourceBuilder = ResourcePoolsBuilder.newResourcePoolsBuilder()
                .heap(100000, EntryUnit.ENTRIES)
                .offheap(300, MemoryUnit.MB);

        cacheManage = CacheManagerBuilder.newCacheManagerBuilder()
                .with(CacheManagerBuilder.persistence(System.getProperty("java.io.tmpdir") + File.separator + ".cache"))
                .withCache(INPUT_SPECTRA_CACHE, CacheConfigurationBuilder
                        .newCacheConfigurationBuilder(String.class, Long.class, spectraResourceBuilder))
                .withCache(CLUSTERS_CACHE, CacheConfigurationBuilder
                        .newCacheConfigurationBuilder(Integer.class, Integer.class, spectraResourceBuilder))
                .withCache(PROTEIN_TO_PSMS_CACHE, CacheConfigurationBuilder
                        .newCacheConfigurationBuilder(String.class, List.class, proteinsResourceBuilder))
                .build();
        cacheManage.init();

        System.out.println("Cache Initialized");
    }

    public static AppCacheManager getInstance() {
        if (instance == null) {
            instance = new AppCacheManager();
        }
        return instance;
    }

    public Cache<Integer, Integer> getClustersCache(){
        return cacheManage.getCache(CLUSTERS_CACHE, Integer.class, Integer.class);
    }

    public Cache<String, Long> getPrideJsonSpectra() {
        return cacheManage.getCache(INPUT_SPECTRA_CACHE, String.class, Long.class);
    }


    public Cache<String, ? extends List> getProteinToPsmsCache() {
        return cacheManage.getCache(PROTEIN_TO_PSMS_CACHE, String.class, List.class);
    }

    public static void closeInstance(){
        if(cacheManage != null)
            cacheManage.close();
        System.out.println("Cache Closed");
    }
}
